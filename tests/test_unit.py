# tests/test_unit.py

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from parser import Node, Material, Section, Element, LoadCase, UniformlyDL, PointLoad, TemperatureL
from element_physics import ElementPhysics
import math_utils

class TestElementPhysics(unittest.TestCase):
    def setUp(self):
        # Setup dummy properties for testing
        self.mat = Material(id="steel", E=2.0e8, alpha=1.2e-5)
        self.sec = Section(id="frame_sec", A=0.01, I=0.0001, d=0.3)
        self.node_i = Node(id=1, x=0.0, y=0.0)
        self.node_j = Node(id=2, x=5.0, y=0.0) 

    def test_local_k_standard_frame(self):
        """Test 1: Compute local [k] for a standard frame element."""
        frame = Element(id="F1", type="frame", node_i=self.node_i, node_j=self.node_j, 
                        material=self.mat, section=self.sec)
        physics = ElementPhysics(frame)
        k_local = physics.get_local_k()

        L = 5.0
        E, I = 2.0e8, 0.0001
        expected_k_v1 = 12 * E * I / (L**3)
        expected_k_r1 = 4 * E * I / L

        self.assertAlmostEqual(k_local[1][1], expected_k_v1, places=3)
        self.assertAlmostEqual(k_local[2][2], expected_k_r1, places=3)

    def test_fef_uniform_distributed_load(self):
        """Test 2: Compute FEFs for a member with a UDL."""
        frame = Element(id="F1", type="frame", node_i=self.node_i, node_j=self.node_j, 
                        material=self.mat, section=self.sec)
        
        udl = UniformlyDL(element=frame, wy=-10.0) 
        fef_local = udl.FEF("fixed-fixed", 5.0)

        w = 10.0  # Use magnitude (abs of -10.0)
        L = 5.0
        
        expected_vy = (w * L) / 2.0
        expected_mz_i = (w * L**2) / 12.0
        
        self.assertAlmostEqual(fef_local[1][0], expected_vy, places=3)
        self.assertAlmostEqual(fef_local[2][0], expected_mz_i, places=3)

    def test_transformation_matrix(self):
        """Test 3: Evaluate coordinate transformation for an inclined element."""
        node_inc_i = Node(id=1, x=0.0, y=0.0)
        node_inc_j = Node(id=2, x=3.0, y=4.0)
        
        frame = Element(id="F2", type="frame", node_i=node_inc_i, node_j=node_inc_j, 
                        material=self.mat, section=self.sec)
        physics = ElementPhysics(frame)
        
        self.assertAlmostEqual(physics.cos_x, 0.6)
        self.assertAlmostEqual(physics.sin_x, 0.8)
        
        fef_local = [[10.0], [0.0], [0.0], [0.0], [0.0], [0.0]]
        k_local_dummy = math_utils.zeros(6, 6)
        
        _, fef_global = physics.transform_to_global(k_local_dummy, fef_local)
        
        self.assertAlmostEqual(fef_global[0][0], 6.0)
        self.assertAlmostEqual(fef_global[1][0], 8.0)
        self.assertAlmostEqual(fef_global[2][0], 0.0)
        
    def test_fef_member_point_load(self):
        """Test 4: Compute FEFs for a point load applied to a member."""
        frame = Element(id="F1", type="frame", node_i=self.node_i, node_j=self.node_j, 
                        material=self.mat, section=self.sec)
        
        P_signed = -20.0
        P = 20.0  # Use magnitude (abs of -20.0)
        a = 2.0
        L = 5.0
        b = L - a
        
        pt_load = PointLoad(element=frame, position=a, fy=P_signed)
        fef_local = pt_load.FEF("fixed-fixed", L)
        
        expected_vy_i = (P * b**2 * (3*a + b)) / (L**3)
        expected_mz_i = (P * a * b**2) / (L**2)
        
        self.assertAlmostEqual(fef_local[1][0], expected_vy_i, places=3)
        self.assertAlmostEqual(fef_local[2][0], expected_mz_i, places=3)

    def test_fef_thermal_gradient(self):
        """Test 5: Compute FEFs for a thermal gradient load (Case 2) - UPDATED per professor's formula.
        
        Originally written with incorrect formula, now updated to reflect:
        delta_T = T_b - T_u (not T_u - T_b)
        M_T = (alpha * delta_T / d) * E * I (exact professor's formula)
        """
        frame = Element(id="F1", type="frame", node_i=self.node_i, node_j=self.node_j, 
                        material=self.mat, section=self.sec)
        
        # Tu = +10, Tb = -10 means delta_T = -20 (bottom colder than top)
        delta_T = 20.0 # Magnitude for reference
        thermal_load = TemperatureL(element=frame, Tu=delta_T/2.0, Tb=-delta_T/2.0)
        fef_local = thermal_load.FEF("fixed-fixed", 5.0)

        # Professor's formula: M_T = (alpha * delta_T / d) * E * I
        # With Tu=10, Tb=-10: delta_T = -10 - 10 = -20
        actual_delta_T = -delta_T  # -20
        expected_M_T = (1.2e-5 * actual_delta_T / 0.3) * (2.0e8 * 0.0001)
        
        # For fixed-fixed: fef[2][0] = -M_T, fef[5][0] = M_T
        self.assertAlmostEqual(fef_local[2][0], -expected_M_T, places=3)
        self.assertAlmostEqual(fef_local[5][0],  expected_M_T, places=3)

    def test_thermal_case1_uniform_truss(self):
        """THERMAL TEST 1: Case 1 - Uniform Temperature Change (Truss).
        
        Verifies: Axial force only for uniformly heated truss.
        Formula: F_T = alpha * T_uniform * E * A
        Expected FEF: [-F_T, 0, 0, F_T, 0, 0]^T (compression when warming)
        """
        # Create a truss element: L=5m, E=2e8, A=0.01, alpha=1.2e-5
        mat_truss = Material(id="steel_truss", E=2.0e8, alpha=1.2e-5)
        sec_truss = Section(id="truss_sec", A=0.01, I=0.0, d=0.0)
        node_i = Node(id=1, x=0.0, y=0.0)
        node_j = Node(id=2, x=5.0, y=0.0)
        
        truss = Element(id="T1", type="truss", node_i=node_i, node_j=node_j,
                       material=mat_truss, section=sec_truss)
        
        # Apply uniform temperature change: dT = 20°C
        thermal_load = TemperatureL(element=truss, Tu=20.0, Tb=20.0)
        fef_local = thermal_load.FEF("pin-pin", 5.0)
        
        # Expected axial force
        alpha = 1.2e-5
        E = 2.0e8
        A = 0.01
        T_uniform = 20.0
        expected_F_T = alpha * T_uniform * E * A
        
        # Verify: Compression force at both ends (sign convention from professor)
        self.assertAlmostEqual(fef_local[0][0], -expected_F_T, places=6)
        self.assertAlmostEqual(fef_local[3][0],  expected_F_T, places=6)
        # Verify: No moment on truss
        self.assertAlmostEqual(fef_local[1][0], 0.0, places=10)
        self.assertAlmostEqual(fef_local[2][0], 0.0, places=10)
        self.assertAlmostEqual(fef_local[4][0], 0.0, places=10)
        self.assertAlmostEqual(fef_local[5][0], 0.0, places=10)

    def test_thermal_case2_gradient_frame(self):
        """THERMAL TEST 2: Case 2 - Temperature Gradient (Frame).
        
        Verifies: Combined axial + bending for gradient-loaded frame.
        Profile: Tu=+10°C (top), Tb=-10°C (bottom), delta_T=20°C
        Formula: M_T = (alpha * delta_T / d) * E * I
        Expected: Axial + moment pair according to professor's decomposition.
        """
        mat_frame = Material(id="concrete", E=3.0e10, alpha=1.0e-5)
        sec_frame = Section(id="beam_sec", A=0.3*0.6, I=0.3*0.6**3/12, d=0.6)
        node_i = Node(id=3, x=0.0, y=0.0)
        node_j = Node(id=4, x=8.0, y=0.0)
        
        frame = Element(id="B1", type="frame", node_i=node_i, node_j=node_j,
                       material=mat_frame, section=sec_frame)
        
        # Gradient load: top warmer than bottom
        thermal_load = TemperatureL(element=frame, Tu=10.0, Tb=-10.0)
        fef_local = thermal_load.FEF("fixed-fixed", 8.0)
        
        # Professor's decomposition
        delta_T = -10.0 - 10.0  # Tb - Tu = -20 (bottom is colder)
        T_uniform = 10.0 + (delta_T / 2.0)  # = 0°C (symmetric)
        
        alpha = 1.0e-5
        E = 3.0e10
        A = 0.3 * 0.6
        I = 0.3 * 0.6**3 / 12
        d = 0.6
        
        F_T = alpha * T_uniform * E * A
        M_T = (alpha * delta_T / d) * E * I
        
        # Verify FEF = [-F_T, 0, -M_T, F_T, 0, M_T]^T
        self.assertAlmostEqual(fef_local[0][0], -F_T, places=4)
        self.assertAlmostEqual(fef_local[1][0],  0.0, places=10)
        self.assertAlmostEqual(fef_local[2][0], -M_T, places=4)
        self.assertAlmostEqual(fef_local[3][0],  F_T, places=4)
        self.assertAlmostEqual(fef_local[4][0],  0.0, places=10)
        self.assertAlmostEqual(fef_local[5][0],  M_T, places=4)

    def test_thermal_case3_combined_frame_pin_fixed(self):
        """THERMAL TEST 3: Case 3 - Combined Action with Release Condition.
        
        Verifies: Moment adjustment when pin release at start (pin-fixed condition).
        Profile: Tu=25°C, Tb=5°C, delta_T=20°C
        Formula: Moment redistributed and shears induced for pin release.
        """
        mat_frame = Material(id="steel_beam", E=2.1e11, alpha=1.2e-5)
        sec_frame = Section(id="box_sec", A=0.04, I=4e-5, d=0.4)
        node_i = Node(id=5, x=0.0, y=0.0)
        node_j = Node(id=6, x=10.0, y=0.0)
        
        # Frame with pin release at start
        frame = Element(id="F3", type="frame", node_i=node_i, node_j=node_j,
                       material=mat_frame, section=sec_frame, release_start=True)
        
        # Combined trapezoidal profile
        thermal_load = TemperatureL(element=frame, Tu=25.0, Tb=5.0)
        fef_local = thermal_load.FEF("pin-fixed", 10.0)
        
        # Professor's decomposition
        delta_T = 5.0 - 25.0  # Tb - Tu = -20 (bottom cooler)
        T_uniform = 25.0 + (delta_T / 2.0)  # = 15°C
        
        alpha = 1.2e-5
        E = 2.1e11
        A = 0.04
        I = 4e-5
        d = 0.4
        L = 10.0
        
        F_T = alpha * T_uniform * E * A
        M_T = (alpha * delta_T / d) * E * I
        
        # For pin-fixed: M_T distributed and shears induced
        # fef[2][0] = 0.0 (pin release)
        # fef[5][0] = M_T - 0.5 * M_T = 0.5 * M_T
        # fef[1][0] = -1.5 * M_T / L, fef[4][0] = 1.5 * M_T / L
        
        self.assertAlmostEqual(fef_local[0][0], -F_T, places=6)
        self.assertAlmostEqual(fef_local[2][0],  0.0, places=10)  # Pin release
        self.assertAlmostEqual(fef_local[3][0],  F_T, places=6)
        self.assertAlmostEqual(fef_local[5][0],  0.5 * M_T, places=4)
        self.assertAlmostEqual(fef_local[1][0], -1.5 * M_T / L, places=4)
        self.assertAlmostEqual(fef_local[4][0],  1.5 * M_T / L, places=4)

if __name__ == '__main__':
    unittest.main()