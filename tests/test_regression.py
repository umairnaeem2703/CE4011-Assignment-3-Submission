# tests/test_regression.py

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from parser import XMLParser
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver
from post_processor import PostProcessor

class TestRegression(unittest.TestCase):

    def setUp(self):
        # Separate tolerances: SAP2000 includes shear deformations by default
        # which causes minor (<1%) deviations from standard Euler-Bernoulli DSM.
        self.disp_tol = 1e-3   # Tolerance for Displacements (m/rad)
        self.force_tol = 0.5   # Tolerance for Forces (kN / kN-m)

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
                if line.startswith("Table:  Element Forces - Frames"):
                    current_table = "forces"
                    continue
                if not line or line.startswith("SAP2000") or line.startswith("Table:"):
                    continue
                    
                tokens = line.split()
                if len(tokens) < 9:
                    continue
                    
                try:
                    # Only process lines where the first column is a numerical ID
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
                    frame_id = tokens[0]
                    station = float(tokens[1])
                    
                    # Tokens: [Frame, Station, Case..., P, V2(Shear), V3, T, M2, M3(Moment)]
                    P = float(tokens[4])
                    V2 = float(tokens[5])
                    M3 = float(tokens[9])
                    
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

        return processor

    # ==========================================================
    # TEST 1 — TRUSS 
    # ==========================================================
    def test_regression_01_truss_displacements_and_reactions(self):
        xml_path = os.path.join(os.path.dirname(__file__), "../data/example1_case1_truss.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/example1_case1_truss_sap2000.txt")
        
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing input or SAP2000 text file")

        sap_disp, sap_react, sap_forces = self._parse_sap2000(sap_path)

        model = XMLParser(xml_path).parse()
        processor = self.run_full_analysis(model)

        disp = processor.displacements
        reactions = processor.reactions
        local_forces = processor.member_forces

        # 1. Check nodal displacement for Support Node 1
        self.assertAlmostEqual(abs(disp[1][0]), abs(sap_disp[1][0]), delta=self.disp_tol)
        self.assertAlmostEqual(abs(disp[1][1]), abs(sap_disp[1][1]), delta=self.disp_tol)

        # 2. Check support reactions for Support Node 1
        self.assertAlmostEqual(abs(reactions[1][0]), abs(sap_react[1][0]), delta=self.force_tol)
        self.assertAlmostEqual(abs(reactions[1][1]), abs(sap_react[1][1]), delta=self.force_tol)
        
        # 3. Check start and end node Axial and Shear for Truss T1
        self.assertIn("T1", local_forces)
        f_T1 = local_forces["T1"]
        sap_f = sap_forces["1"]
        
        # Truss local vector: [Fx_i, Fy_i, Fx_j, Fy_j]^T
        self.assertAlmostEqual(abs(f_T1[0][0]), abs(sap_f['i'][0]), delta=self.force_tol) # Axial i
        self.assertAlmostEqual(abs(f_T1[2][0]), abs(sap_f['j'][0]), delta=self.force_tol) # Axial j
        
        self.assertAlmostEqual(abs(f_T1[1][0]), abs(sap_f['i'][1]), delta=self.force_tol) # Shear i (0)
        self.assertAlmostEqual(abs(f_T1[3][0]), abs(sap_f['j'][1]), delta=self.force_tol) # Shear j (0)

    # ==========================================================
    # TEST 2 — FRAME 
    # ==========================================================
    def test_regression_02_frame_full_validation(self):
        xml_path = os.path.join(os.path.dirname(__file__), "../data/example2_frame.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/example2_frame_sap2000.txt")
        
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing input or SAP2000 text file")

        sap_disp, sap_react, sap_forces = self._parse_sap2000(sap_path)

        model = XMLParser(xml_path).parse()
        processor = self.run_full_analysis(model)

        disp = processor.displacements
        reactions = processor.reactions
        local_forces = processor.member_forces

        # 1. Check nodal displacement for Support Node 1
        self.assertAlmostEqual(abs(disp[1][0]), abs(sap_disp[1][0]), delta=self.disp_tol)
        self.assertAlmostEqual(abs(disp[1][1]), abs(sap_disp[1][1]), delta=self.disp_tol)
        self.assertAlmostEqual(abs(disp[1][2]), abs(sap_disp[1][2]), delta=self.disp_tol)

        # 2. Check support reactions for Support Node 1
        self.assertAlmostEqual(abs(reactions[1][0]), abs(sap_react[1][0]), delta=self.force_tol)
        self.assertAlmostEqual(abs(reactions[1][1]), abs(sap_react[1][1]), delta=self.force_tol)
        self.assertAlmostEqual(abs(reactions[1][2]), abs(sap_react[1][2]), delta=self.force_tol)

        # 3. Check internal forces for Frame F1
        self.assertIn("F1", local_forces)
        f_F1 = local_forces["F1"]
        sap_f = sap_forces["1"]
        
        # Frame local vector: [Fx_i, Vy_i, Mz_i, Fx_j, Vy_j, Mz_j]^T
        self.assertAlmostEqual(abs(f_F1[0][0]), abs(sap_f['i'][0]), delta=self.force_tol) # Axial i
        self.assertAlmostEqual(abs(f_F1[3][0]), abs(sap_f['j'][0]), delta=self.force_tol) # Axial j
        
        self.assertAlmostEqual(abs(f_F1[1][0]), abs(sap_f['i'][1]), delta=self.force_tol) # Shear i
        self.assertAlmostEqual(abs(f_F1[4][0]), abs(sap_f['j'][1]), delta=self.force_tol) # Shear j
        
        self.assertAlmostEqual(abs(f_F1[2][0]), abs(sap_f['i'][2]), delta=self.force_tol) # Moment i
        self.assertAlmostEqual(abs(f_F1[5][0]), abs(sap_f['j'][2]), delta=self.force_tol) # Moment j

    # ==========================================================
    # TEST 3 — MIXED FRAME-TRUSS 
    # ==========================================================
    def test_regression_03_mixed_structure_results(self):
        xml_path = os.path.join(os.path.dirname(__file__), "../data/example3_frame_truss.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/example3_frame_truss_sap2000.txt")
        
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing input or SAP2000 text file")

        sap_disp, sap_react, sap_forces = self._parse_sap2000(sap_path)

        model = XMLParser(xml_path).parse()
        processor = self.run_full_analysis(model)

        disp = processor.displacements
        reactions = processor.reactions
        local_forces = processor.member_forces

        # 1. Check nodal displacement for Support Node 1
        self.assertAlmostEqual(abs(disp[1][0]), abs(sap_disp[1][0]), delta=self.disp_tol)
        self.assertAlmostEqual(abs(disp[1][1]), abs(sap_disp[1][1]), delta=self.disp_tol)
        self.assertAlmostEqual(abs(disp[1][2]), abs(sap_disp[1][2]), delta=self.disp_tol)

        # 2. Check support reactions for Support Node 1
        self.assertAlmostEqual(abs(reactions[1][0]), abs(sap_react[1][0]), delta=self.force_tol)
        self.assertAlmostEqual(abs(reactions[1][1]), abs(sap_react[1][1]), delta=self.force_tol)
        self.assertAlmostEqual(abs(reactions[1][2]), abs(sap_react[1][2]), delta=self.force_tol)

        # 3. Check internal forces for Frame F1
        self.assertIn("F1", local_forces)
        f_F1 = local_forces["F1"]
        sap_f = sap_forces["1"]
        
        self.assertAlmostEqual(abs(f_F1[0][0]), abs(sap_f['i'][0]), delta=self.force_tol) # Axial i
        self.assertAlmostEqual(abs(f_F1[3][0]), abs(sap_f['j'][0]), delta=self.force_tol) # Axial j
        
        self.assertAlmostEqual(abs(f_F1[1][0]), abs(sap_f['i'][1]), delta=self.force_tol) # Shear i
        self.assertAlmostEqual(abs(f_F1[4][0]), abs(sap_f['j'][1]), delta=self.force_tol) # Shear j
        
        self.assertAlmostEqual(abs(f_F1[2][0]), abs(sap_f['i'][2]), delta=self.force_tol) # Moment i
        self.assertAlmostEqual(abs(f_F1[5][0]), abs(sap_f['j'][2]), delta=self.force_tol) # Moment j

if __name__ == '__main__':
    unittest.main()