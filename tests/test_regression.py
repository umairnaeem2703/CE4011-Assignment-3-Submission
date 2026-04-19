# tests/test_regression.py

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from parser import XMLParser, StructuralModel, Node, Element, Material, Section, Support, LoadCase, NodalLoad
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver
from post_processor import PostProcessor

class TestRegression(unittest.TestCase):

    def setUp(self):
        self.disp_tol = 1e-3   # Standard Tolerance for SAP2000 Displacements
        self.force_tol = 0.5   # Standard Tolerance for SAP2000 Forces
        self.exact_tol = 1e-5  # Stricter Tolerance for Analytical/Hardcoded checks

    def _parse_sap2000(self, xml_path):
        """Helper method to dynamically parse SAP2000 results from a text file."""
        disp = {}
        react = {}
        forces = {}
        current_table = None
        
        with open(xml_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("Table:  Joint Displacements"):
                    current_table = "disp"
                    continue
                if line.startswith("Table:  Joint Reactions"):
                    current_table = "react"
                    continue
                if line.startswith("Table:  Element Joint Forces - Frames"):
                    current_table = "forces"
                    continue
                if not line or line.startswith("SAP2000") or line.startswith("Table:"):
                    continue
                    
                tokens = line.split()
                if len(tokens) < 9:
                    continue
                
                # For disp and react, verify first token is numeric
                if current_table in ("disp", "react"):
                    try:
                        first_val = float(tokens[0])
                    except ValueError:
                        continue
                    
                if current_table == "disp":
                    node_id = int(tokens[0])
                    # SAP2000 (XZ Plane): U1=X, U3=Y, R2=Rot_Y (Opposite to our RZ)
                    disp[node_id] = (float(tokens[3]), float(tokens[5]), -float(tokens[7]))
                    
                elif current_table == "react":
                    node_id = int(tokens[0])
                    # SAP2000 (XZ Plane): F1=Fx, F3=Fy, M2=Mz (Opposite to our Mz)
                    react[node_id] = (float(tokens[3]), float(tokens[5]), -float(tokens[7]))
                    
                elif current_table == "forces" and len(tokens) >= 10:
                    try:
                        frame_id = tokens[0]
                        station = float(tokens[1])
                        
                        # SAP2000 XZ-plane coordinate system: X=1(horizontal), Z=3(vertical), Y=2(out-of-plane)
                        # For 2D analysis in our XY plane (with Y as vertical), map: SAP X->our X, SAP Z->our Y
                        P = float(tokens[4])    # F1 - Axial (X direction)
                        V2 = float(tokens[6])   # F3 - Vertical shear (Z direction, our F3)
                        M3 = -float(tokens[8])   # M2 - In-plane moment (negated: SAP convention opposite to ours)
                    except (ValueError, IndexError):
                        continue
                    
                    if frame_id not in forces:
                        forces[frame_id] = {'min_st': station, 'max_st': station, 'i': (P, V2, M3), 'j': (P, V2, M3)}
                    else:
                        if station < forces[frame_id]['min_st']:
                            forces[frame_id]['min_st'] = station
                            forces[frame_id]['i'] = (P, V2, M3)
                        if station > forces[frame_id]['max_st']:
                            forces[frame_id]['max_st'] = station
                            forces[frame_id]['j'] = (P, V2, M3)
                            
        return disp, react, forces

    def run_full_analysis(self, model, load_case="LC1"):
        """Reusable pipeline"""
        opt = DOFOptimizer(model)
        num_eq, semi_bw, _ = opt.optimize()

        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K_banded, F_global = assembler.assemble(load_case)

        solver = BandedSolver(K_banded, F_global, semi_bw)
        D_active = solver.solve()

        processor = PostProcessor(model, D_active, load_case)

        return processor, F_global

    # ==========================================================
    # TEST 1 — TRUSS (SAP2000 Benchmark)
    # ==========================================================
    def test_regression_01_truss_displacements_and_reactions(self):
        xml_path = os.path.join(os.path.dirname(__file__), "../data/example1_case1_truss.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/example1_case1_truss_sap2000.txt")
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing input or SAP2000 text file")

        sap_disp, sap_react, sap_forces = self._parse_sap2000(sap_path)
        model = XMLParser(xml_path).parse()
        processor, _ = self.run_full_analysis(model)

        disp = processor.displacements
        reactions = processor.reactions
        local_forces = processor.member_forces

        self.assertAlmostEqual(disp[1][0], sap_disp[1][0], delta=self.disp_tol)
        self.assertAlmostEqual(disp[1][1], sap_disp[1][1], delta=self.disp_tol)
        self.assertAlmostEqual(reactions[1][0], sap_react[1][0], delta=self.force_tol)
        self.assertAlmostEqual(reactions[1][1], sap_react[1][1], delta=self.force_tol)
        self.assertIn("T1", local_forces)
        
        # SAP2000 reports GLOBAL forces; convert to local for comparison
        el_T1 = model.elements["T1"]
        dx = el_T1.node_j.x - el_T1.node_i.x
        dy = el_T1.node_j.y - el_T1.node_i.y
        L = (dx**2 + dy**2)**0.5
        cos_x, sin_x = dx / L, dy / L
        
        # Transform SAP2000 global forces to local coordinates
        sap_f_global = sap_forces["T1"]
        sap_P_i = sap_f_global['i'][0] * cos_x + sap_f_global['i'][1] * sin_x
        sap_P_j = sap_f_global['j'][0] * cos_x + sap_f_global['j'][1] * sin_x
        
        # Compare local forces
        f_T1 = local_forces["T1"]
        self.assertAlmostEqual(f_T1[0][0], sap_P_i, delta=self.force_tol)
        self.assertAlmostEqual(f_T1[2][0], sap_P_j, delta=self.force_tol)

    # ==========================================================
    # TEST 2 — FRAME (SAP2000 Benchmark)
    # ==========================================================
    def test_regression_02_frame_full_validation(self):
        xml_path = os.path.join(os.path.dirname(__file__), "../data/example2_frame.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/example2_frame_sap2000.txt")
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing input or SAP2000 text file")

        sap_disp, sap_react, sap_forces = self._parse_sap2000(sap_path)
        model = XMLParser(xml_path).parse()
        processor, _ = self.run_full_analysis(model)

        disp = processor.displacements
        reactions = processor.reactions
        local_forces = processor.member_forces

        # Validate displacements
        self.assertAlmostEqual(disp[1][0], sap_disp[1][0], delta=self.disp_tol)
        self.assertAlmostEqual(disp[1][1], sap_disp[1][1], delta=self.disp_tol)
        self.assertAlmostEqual(disp[1][2], sap_disp[1][2], delta=self.disp_tol)
        
        # Validate reactions
        self.assertAlmostEqual(reactions[1][0], sap_react[1][0], delta=self.force_tol)
        self.assertAlmostEqual(reactions[1][1], sap_react[1][1], delta=self.force_tol)
        self.assertAlmostEqual(reactions[1][2], sap_react[1][2], delta=self.force_tol)

        # Validate member force structure (6 components: P, V, M at i and j ends)
        self.assertIn("F1", local_forces)
        f_F1 = local_forces["F1"]
        self.assertEqual(len(f_F1), 6, "Frame element should have 6 force components")
        
        # Verify forces are computed (non-trivial values)
        all_forces = [f[0] for f in f_F1]
        self.assertTrue(any(abs(f) > 0.1 for f in all_forces), "At least some forces should be significant")

    # ==========================================================
    # TEST 3 — MIXED FRAME-TRUSS (SAP2000 Benchmark)
    # ==========================================================
    def test_regression_03_mixed_structure_results(self):
        xml_path = os.path.join(os.path.dirname(__file__), "../data/example3_frame_truss.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/example3_frame_truss_sap2000.txt")
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing input or SAP2000 text file")

        sap_disp, sap_react, sap_forces = self._parse_sap2000(sap_path)
        model = XMLParser(xml_path).parse()
        processor, _ = self.run_full_analysis(model)

        disp = processor.displacements
        reactions = processor.reactions
        local_forces = processor.member_forces

        # CRITICAL FIX: Check SIGNED values to verify load directionality
        # A downward load (fy=-20) must produce downward displacement (negative UY)
        self.assertAlmostEqual(disp[1][0], sap_disp[1][0], delta=self.disp_tol)
        self.assertAlmostEqual(disp[1][1], sap_disp[1][1], delta=self.disp_tol)
        self.assertAlmostEqual(disp[1][2], sap_disp[1][2], delta=self.disp_tol)
        self.assertAlmostEqual(reactions[1][0], sap_react[1][0], delta=self.force_tol)
        self.assertAlmostEqual(reactions[1][1], sap_react[1][1], delta=self.force_tol)
        self.assertAlmostEqual(reactions[1][2], sap_react[1][2], delta=self.force_tol)

        self.assertIn("F1", local_forces)
        f_F1 = local_forces["F1"]
        sap_f = sap_forces["F1"]
        self.assertAlmostEqual(f_F1[0][0], sap_f['i'][0], delta=self.force_tol)
        self.assertAlmostEqual(f_F1[3][0], sap_f['j'][0], delta=self.force_tol)
        self.assertAlmostEqual(f_F1[1][0], sap_f['i'][1], delta=self.force_tol)
        self.assertAlmostEqual(f_F1[4][0], sap_f['j'][1], delta=self.force_tol)
        self.assertAlmostEqual(f_F1[2][0], sap_f['i'][2], delta=self.force_tol)
        self.assertAlmostEqual(f_F1[5][0], sap_f['j'][2], delta=self.force_tol)

    # ==========================================================
    # TEST 4 — EG1 TRUSS TEMPERATURE (Hardcoded Answers)
    # ==========================================================
    def test_regression_04_eg1_truss_temp(self):
        """Validates the 3-truss temperature problem against exact whiteboard answers."""
        xml_path = os.path.join(os.path.dirname(__file__), "../data/eg1_truss_temp.xml")
        if not os.path.exists(xml_path):
            self.skipTest("Missing input file eg1_truss_temp.xml")

        model = XMLParser(xml_path).parse()
        processor, F_global = self.run_full_analysis(model)

        # Expected Vector F (Equivalent Nodal Loads derived from FEF)
        # N1 is free in X and Y.
        n1 = model.nodes[1]
        self.assertAlmostEqual(F_global[n1.dofs[0]][0], -72.08, delta=0.01)
        self.assertAlmostEqual(F_global[n1.dofs[1]][0], -32.23, delta=0.01)

        # Expected Displacements u (in metres)
        disp = processor.displacements
        self.assertAlmostEqual(disp[1][0], -0.000405, delta=self.exact_tol)
        self.assertAlmostEqual(disp[1][1], -0.000070, delta=self.exact_tol)

    # ==========================================================
    # TEST 5 — EG2 BEAM TEMPERATURE (Hardcoded Answers)
    # ==========================================================
    def test_regression_05_eg2_beam_temp(self):
        """Validates the 2-span beam thermal gradient problem against exact whiteboard answers."""
        xml_path = os.path.join(os.path.dirname(__file__), "../data/eg2_beam_temp.xml")
        if not os.path.exists(xml_path):
            self.skipTest("Missing input file eg2_beam_temp.xml")

        model = XMLParser(xml_path).parse()
        processor, F_global = self.run_full_analysis(model)

        # Expected Vector F (Equivalent Moments at nodes 1 and 2)
        n1 = model.nodes[1]
        n2 = model.nodes[2]
        self.assertAlmostEqual(F_global[n1.dofs[2]][0],  12.0, delta=0.01)
        self.assertAlmostEqual(F_global[n2.dofs[2]][0],  -6.0, delta=0.01)

        # Expected Displacements u (Rotations in radians)
        disp = processor.displacements
        self.assertAlmostEqual(disp[1][2],  0.00086, delta=self.exact_tol)
        self.assertAlmostEqual(disp[2][2], -0.00052, delta=self.exact_tol)

if __name__ == '__main__':
    unittest.main()