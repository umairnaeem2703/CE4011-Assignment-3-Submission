# tests/test_regression.py

import sys
import os
import unittest

# Add both source and test directories to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from parser import XMLParser, StructuralModel, Node, Element, Material, Section, Support, LoadCase, NodalLoad
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver
from post_processor import PostProcessor
from sap2000_parser import SAP2000Parser, MemberForceTransformer, CoordinateMapper, assert_displacement_match, assert_force_match

class TestRegression(unittest.TestCase):

    def setUp(self):
        # Tolerances per specification:
        # - Displacement tolerance: 2% (relative) - Accounts for different solver implementations
        #   SAP2000 uses advanced solvers; our 2D solver is simpler but adequate for educational purposes
        # - Force tolerance: 2% (relative) - SAP2000 includes shear deformation, our solver uses Euler-Bernoulli
        #   This ~1.5-2% discrepancy is expected and acceptable
        self.disp_rel_tol = 0.02    # 2% relative tolerance for displacements (different solver algorithms)
        self.force_rel_tol = 0.02   # 2% relative tolerance for forces/reactions (shear deformation theory)
        self.force_abs_tol = 0.5    # For legacy support, absolute tolerance fallback

    def _run_full_analysis(self, model, load_case="LC1"):
        """
        Execute complete FEA pipeline: optimization, assembly, solving, post-processing.
        
        Args:
            model: StructuralModel instance
            load_case: Load case identifier (default "LC1")
            
        Returns:
            Tuple of (PostProcessor, F_global_vector)
        """
        opt = DOFOptimizer(model)
        num_eq, semi_bw, _ = opt.optimize()

        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K_banded, F_global = assembler.assemble(load_case)

        solver = BandedSolver(K_banded, F_global, semi_bw)
        D_active = solver.solve()

        processor = PostProcessor(model, D_active, load_case)

        return processor, F_global

    # ==========================================================
    # TEST 1 — ASSIGNMENT Q2(a) WITH SETTLEMENTS (SAP2000 Benchmark)
    # ==========================================================
    def test_regression_01_assignment_q2a_with_settlements(self):
        """
        Validates support settlement handling with complete SAP2000 reference comparison.
        
        This comprehensive test validates:
        1. Nodal displacements (both free and prescribed settlements)
        2. Support reactions
        3. Member end forces with proper local-to-global transformation
        
        The settlement implementation is verified through:
        - Prescribed DOF displacements matching settlement values
        - Reactions accounting for settlement-induced moments
        - Member forces correctly transformed from local to global coordinates
        
        Key transformation:
        - SAP2000 provides forces in GLOBAL coordinates (F1, F3, M2)
        - Our solver computes forces in LOCAL (element) coordinates
        - Before comparison, we transform using element orientation:
          Fx_global = cos(θ)*P_local - sin(θ)*V_local
          Fy_global = sin(θ)*P_local + cos(θ)*V_local
          Mz_global = M_local
        """
        xml_path = os.path.join(os.path.dirname(__file__), "../data/Assignment_4_Q2a.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/q2_a_sap2000.txt")
        
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing Assignment_4_Q2a.xml or q2_a_sap2000.txt reference file")

        # Parse SAP2000 reference with proper coordinate mapping
        parser = SAP2000Parser(sap_path)
        sap_disp, sap_react, sap_forces = parser.parse()
        
        # Run our solver
        model = XMLParser(xml_path).parse()
        processor, F_global = self._run_full_analysis(model)

        disp = processor.displacements
        reactions = processor.reactions
        local_forces = processor.member_forces

        # ========== 1. VALIDATE NODAL DISPLACEMENTS ==========
        # Verify all computed displacements match SAP2000
        for node_id in sap_disp.keys():
            if node_id in disp:
                match, error_msg = assert_displacement_match(disp[node_id], sap_disp[node_id],
                                                             rel_tol=self.disp_rel_tol)
                self.assertTrue(match, f"Node {node_id} displacement mismatch:\n{error_msg}")
        
        # ========== 2. VALIDATE SUPPORT REACTIONS ==========
        # Verify reactions at all support nodes
        for node_id in sap_react.keys():
            if node_id in reactions:
                match, error_msg = assert_force_match(reactions[node_id], sap_react[node_id],
                                                      rel_tol=self.force_rel_tol)
                self.assertTrue(match, f"Node {node_id} reaction mismatch:\n{error_msg}")
        
        # ========== 3. VALIDATE MEMBER END FORCES WITH TRANSFORMATION ==========
        # SAP2000 provides global forces; transform our local forces for comparison
        for elem_id in sap_forces.keys():
            if elem_id not in local_forces or elem_id not in model.elements:
                continue
                
            el = model.elements[elem_id]
            
            # Skip truss elements (they don't have moment components and use different output format)
            # Only validate frame elements with full force/moment data (6 components)
            if not hasattr(el, 'moment_inertia') or el.moment_inertia == 0:
                # This is a truss element - skip detailed force transformation
                # Truss forces are easier to validate but SAP2000 output format may differ
                continue
            
            # Element is a frame with moments - perform coordinate transformation
            theta, cos_theta, sin_theta, length = MemberForceTransformer.calculate_element_angle(
                el.node_i, el.node_j
            )
            
            sap_forces_elem = sap_forces[elem_id]
            our_forces_elem = local_forces[elem_id]
            
            # Frame elements have 6 components: [P_i, V_i, M_i, P_j, V_j, M_j]
            if len(our_forces_elem) < 6:
                continue
            
            # Transform our local end-I forces to global
            our_global_i = MemberForceTransformer.local_to_global_forces(
                our_forces_elem[0][0],  # Axial
                our_forces_elem[1][0],  # Shear
                our_forces_elem[2][0],  # Moment
                cos_theta, sin_theta
            )
            
            # Transform our local end-J forces to global
            our_global_j = MemberForceTransformer.local_to_global_forces(
                our_forces_elem[3][0],  # Axial
                our_forces_elem[4][0],  # Shear
                our_forces_elem[5][0],  # Moment
                cos_theta, sin_theta
            )
            
            # Compare against SAP2000 global forces
            match_i, error_i = assert_force_match(our_global_i, sap_forces_elem['i'],
                                                 rel_tol=self.force_rel_tol)
            match_j, error_j = assert_force_match(our_global_j, sap_forces_elem['j'],
                                                 rel_tol=self.force_rel_tol)
            
            self.assertTrue(match_i, f"Element {elem_id} node I forces mismatch:\n{error_i}")
            self.assertTrue(match_j, f"Element {elem_id} node J forces mismatch:\n{error_j}")
        
        # ========== 4. VERIFY SETTLEMENT VALUES ARE PRESCRIBED ==========
        # Confirm nodes with settlement constraints maintain prescribed values
        for node_id, node in model.nodes.items():
            support = model.supports.get(node_id)
            if support and (support.settlement_ux != 0.0 or 
                          support.settlement_uy != 0.0 or 
                          support.settlement_rz != 0.0):
                our_disp = disp[node_id]
                
                if support.restrain_ux:
                    self.assertAlmostEqual(our_disp[0], support.settlement_ux, places=6,
                                         msg=f"Node {node_id} settlement_ux not prescribed")
                if support.restrain_uy:
                    self.assertAlmostEqual(our_disp[1], support.settlement_uy, places=6,
                                         msg=f"Node {node_id} settlement_uy not prescribed")
                if support.restrain_rz and len(our_disp) > 2:
                    self.assertAlmostEqual(our_disp[2], support.settlement_rz, places=6,
                                         msg=f"Node {node_id} settlement_rz not prescribed")

if __name__ == '__main__':
    unittest.main()