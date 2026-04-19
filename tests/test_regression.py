# tests/test_regression.py

import sys
import os
import re
import unittest
import math

# Add both source and test directories to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from parser import XMLParser
from sap2000_parser import SAP2000Parser, MemberForceTransformer, assert_displacement_match, assert_force_match


class SolverResultsParser:
    """
    Parse solver-generated report text files containing:
    - NODAL DISPLACEMENTS
    - MEMBER LOCAL END FORCES
    - SUPPORT REACTIONS
    """
    
    @staticmethod
    def parse(report_path):
        """
        Parse solver results report.
        
        Returns:
            Tuple of (displacements, reactions, local_forces)
            - displacements: {node_id: (ux, uy, rz)}
            - reactions: {node_id: (fx, fy, mz)}
            - local_forces: {elem_id: {'i': (axial, shear, moment), 'j': (axial, shear, moment)}}
        """
        with open(report_path, "r") as f:
            lines = f.readlines()
        
        displacements = {}
        reactions = {}
        local_forces_by_node = {}
        
        section = None
        current_elem = None
        
        for line in lines:
            stripped = line.strip()
            
            if not stripped or stripped.startswith("=" * 10) or stripped.startswith("-" * 10):
                continue
            
            if "NODAL DISPLACEMENTS" in stripped:
                section = "disp"
                continue
            elif "MEMBER LOCAL END FORCES" in stripped:
                section = "forces"
                continue
            elif "SUPPORT REACTIONS" in stripped:
                section = "react"
                continue
            
            if section == "disp" and "Node" not in stripped:
                parts = stripped.split()
                if len(parts) >= 4:
                    try:
                        node_id = int(parts[0])
                        ux = float(parts[1])
                        uy = float(parts[2])
                        rz = float(parts[3])
                        displacements[node_id] = (ux, uy, rz)
                    except (ValueError, IndexError):
                        pass
            
            elif section == "forces" and "Element" not in stripped and "Node" not in stripped:
                parts = stripped.split()
                if len(parts) < 4:
                    continue
                
                try:
                    if parts[0] and not parts[0][0].isdigit():
                        current_elem = parts[0]
                        node_id = int(parts[1])
                        axial = float(parts[2])
                        shear = float(parts[3])
                        moment = float(parts[4])
                    elif current_elem and parts[0].isdigit():
                        node_id = int(parts[0])
                        axial = float(parts[1])
                        shear = float(parts[2])
                        moment = float(parts[3])
                    else:
                        continue
                    
                    if current_elem not in local_forces_by_node:
                        local_forces_by_node[current_elem] = {}
                    local_forces_by_node[current_elem][node_id] = (axial, shear, moment)
                except (ValueError, IndexError):
                    pass
            
            elif section == "react" and "Node" not in stripped:
                parts = stripped.split()
                if len(parts) >= 3:
                    try:
                        node_id = int(parts[0])
                        fx = float(parts[1])
                        fy = float(parts[2])
                        mz = 0.0 if (len(parts) < 4 or "Free" in parts[3] or "free" in parts[3]) else float(parts[3])
                        reactions[node_id] = (fx, fy, mz)
                    except (ValueError, IndexError):
                        pass
        
        local_forces = {}
        for elem_id, node_map in local_forces_by_node.items():
            if len(node_map) >= 2:
                node_ids = sorted(node_map.keys())
                local_forces[elem_id] = {
                    'i': node_map[node_ids[0]],
                    'j': node_map[node_ids[1]]
                }
        
        return displacements, reactions, local_forces


class TestRegression(unittest.TestCase):

    def setUp(self):
        # Tolerances:
        # - Displacement: 2% (different solver implementations)
        # - Forces/Reactions: 2% (Euler-Bernoulli vs other formulations)
        self.disp_rel_tol = 0.02
        self.force_rel_tol = 0.02

    def _get_element_angle(self, node_i_id, node_j_id, model):
        """Get angle and direction cosines for an element."""
        node_i = model.nodes[node_i_id]
        node_j = model.nodes[node_j_id]
        dx = node_j.x - node_i.x
        dy = node_j.y - node_i.y
        length = (dx**2 + dy**2)**0.5
        cos_theta = dx / length if length > 0 else 1.0
        sin_theta = dy / length if length > 0 else 0.0
        return cos_theta, sin_theta

    def _compare_results(self, sap_disp, sap_react, sap_forces, solver_disp, solver_react, solver_local_forces, model, test_name):
        """
        Compare SAP2000 results with solver results.
        
        For member forces:
        - SAP2000 provides global forces (F1, F3, M2)
        - Solver provides local forces (axial, shear, moment)
        - Transform solver local to global for comparison
        - For inclined members: compare F1/F3 resultant with transformed axial/shear
        """
        # VALIDATE DISPLACEMENTS
        for node_id in sap_disp.keys():
            if node_id in solver_disp:
                match, error_msg = assert_displacement_match(
                    solver_disp[node_id],
                    sap_disp[node_id],
                    rel_tol=self.disp_rel_tol
                )
                self.assertTrue(match, f"{test_name} node {node_id} displacement mismatch:\n{error_msg}")

        # VALIDATE REACTIONS
        for node_id in sap_react.keys():
            if node_id in solver_react:
                match, error_msg = assert_force_match(
                    solver_react[node_id],
                    sap_react[node_id],
                    rel_tol=self.force_rel_tol
                )
                self.assertTrue(match, f"{test_name} node {node_id} reaction mismatch:\n{error_msg}")

        # VALIDATE MEMBER END FORCES WITH LOCAL-TO-GLOBAL TRANSFORMATION
        for elem_id in sap_forces.keys():
            if elem_id not in solver_local_forces or elem_id not in model.elements:
                continue

            el = model.elements[elem_id]
            cos_theta, sin_theta = self._get_element_angle(el.node_i.id, el.node_j.id, model)

            solver_forces_elem = solver_local_forces[elem_id]
            sap_forces_elem = sap_forces[elem_id]

            # Transform solver local forces to global
            solver_global_i = MemberForceTransformer.local_to_global_forces(
                solver_forces_elem['i'][0],  # axial
                solver_forces_elem['i'][1],  # shear
                solver_forces_elem['i'][2],  # moment
                cos_theta,
                sin_theta,
            )

            solver_global_j = MemberForceTransformer.local_to_global_forces(
                solver_forces_elem['j'][0],  # axial
                solver_forces_elem['j'][1],  # shear
                solver_forces_elem['j'][2],  # moment
                cos_theta,
                sin_theta,
            )

            match_i, error_i = assert_force_match(
                solver_global_i,
                sap_forces_elem['i'],
                rel_tol=self.force_rel_tol,
            )
            match_j, error_j = assert_force_match(
                solver_global_j,
                sap_forces_elem['j'],
                rel_tol=self.force_rel_tol,
            )

            self.assertTrue(match_i, f"{test_name} element {elem_id} node I forces mismatch:\n{error_i}")
            self.assertTrue(match_j, f"{test_name} element {elem_id} node J forces mismatch:\n{error_j}")

    # ==========================================================
    # TEST 1 — ASSIGNMENT Q2(a) WITH SETTLEMENTS
    # ==========================================================
    def test_regression_01_assignment_q2a(self):
        """
        Validate Q2a solver results against SAP2000 reference.
        
        Compares:
        - Nodal displacements
        - Support reactions
        - Member end forces (with local-to-global transformation for inclined members)
        """
        xml_path = os.path.join(os.path.dirname(__file__), "../data/Assignment_4_Q2a.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/q2_a_sap2000.txt")
        solver_path = os.path.join(os.path.dirname(__file__), "../results/Assignment_4_Q2a_LC1_results.txt")

        if not os.path.exists(sap_path) or not os.path.exists(solver_path):
            self.skipTest("Missing q2_a_sap2000.txt or Assignment_4_Q2a_LC1_results.txt")

        model = XMLParser(xml_path).parse()

        sap_parser = SAP2000Parser(sap_path)
        sap_disp, sap_react, sap_forces = sap_parser.parse()

        solver_disp, solver_react, solver_local_forces = SolverResultsParser.parse(solver_path)

        self._compare_results(sap_disp, sap_react, sap_forces, 
                            solver_disp, solver_react, solver_local_forces,
                            model, "Q2a")

    # ==========================================================
    # TEST 2 — ASSIGNMENT Q2(b) THERMAL LOADING
    # ==========================================================
    def test_regression_02_assignment_q2b(self):
        """
        Validate Q2b solver results against SAP2000 reference.
        
        Compares:
        - Nodal displacements
        - Support reactions
        - Member end forces (with local-to-global transformation for inclined members)
        """
        xml_path = os.path.join(os.path.dirname(__file__), "../data/Assignment_4_Q2b.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/q2_b_sap2000.txt")
        solver_path = os.path.join(os.path.dirname(__file__), "../results/Assignment_4_Q2b_LC1_results.txt")

        if not os.path.exists(sap_path) or not os.path.exists(solver_path):
            self.skipTest("Missing q2_b_sap2000.txt or Assignment_4_Q2b_LC1_results.txt")

        model = XMLParser(xml_path).parse()

        sap_parser = SAP2000Parser(sap_path)
        sap_disp, sap_react, sap_forces = sap_parser.parse()

        solver_disp, solver_react, solver_local_forces = SolverResultsParser.parse(solver_path)

        self._compare_results(sap_disp, sap_react, sap_forces,
                            solver_disp, solver_react, solver_local_forces,
                            model, "Q2b")

if __name__ == '__main__':
    unittest.main()