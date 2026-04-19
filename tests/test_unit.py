# tests/test_unit.py

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from parser import Node, Material, Section, Element, Support, TemperatureL, UniformlyDL
from element_physics import ElementPhysics
import math_utils


class TestElementPhysics(unittest.TestCase):
    """Core element physics tests."""
    
    def setUp(self):
        self.mat = Material(id="steel", E=2.0e8, alpha=1.2e-5)
        self.sec = Section(id="frame_sec", A=0.01, I=0.0001, d=0.3)
        self.node_i = Node(id=1, x=0.0, y=0.0)
        self.node_j = Node(id=2, x=5.0, y=0.0)

    def test_local_stiffness_matrix_frame(self):
        """Verify local stiffness matrix computation for frame element."""
        frame = Element(id="F1", type="frame", node_i=self.node_i, node_j=self.node_j,
                       material=self.mat, section=self.sec)
        physics = ElementPhysics(frame)
        k_local = physics.get_local_k()

        L = 5.0
        E, I = 2.0e8, 0.0001
        expected_k_v = 12 * E * I / (L**3)
        expected_k_r = 4 * E * I / L

        self.assertAlmostEqual(k_local[1][1], expected_k_v, places=3)
        self.assertAlmostEqual(k_local[2][2], expected_k_r, places=3)

    def test_element_transformation_inclined(self):
        """Verify coordinate transformation for inclined element (3-4-5 triangle)."""
        node_i = Node(id=1, x=0.0, y=0.0)
        node_j = Node(id=2, x=3.0, y=4.0)
        
        frame = Element(id="F2", type="frame", node_i=node_i, node_j=node_j,
                       material=self.mat, section=self.sec)
        physics = ElementPhysics(frame)
        
        self.assertAlmostEqual(physics.cos_x, 0.6, places=6)
        self.assertAlmostEqual(physics.sin_x, 0.8, places=6)
        
        fef_local = [[10.0], [0.0], [0.0], [0.0], [0.0], [0.0]]
        k_local_dummy = math_utils.zeros(6, 6)
        
        _, fef_global = physics.transform_to_global(k_local_dummy, fef_local)
        
        self.assertAlmostEqual(fef_global[0][0], 6.0, places=6)
        self.assertAlmostEqual(fef_global[1][0], 8.0, places=6)


class TestThermalLoading(unittest.TestCase):
    """Unit tests for thermal loading fixed-end forces."""
    
    def setUp(self):
        self.mat_frame = Material(id="concrete", E=3.0e10, alpha=1.0e-5)
        self.sec_frame = Section(id="beam_sec", A=0.18, I=0.0054, d=0.6)
        self.node_i = Node(id=1, x=0.0, y=0.0)
        self.node_j = Node(id=2, x=8.0, y=0.0)

    def test_thermal_uniform_temperature_truss(self):
        """TEST 1: Uniform temperature change in truss.
        
        Verifies axial thermal force only (no bending).
        Formula: F_T = alpha * T_uniform * E * A
        Expected: [-F_T, 0, 0, F_T, 0, 0]^T
        """
        mat_truss = Material(id="steel", E=2.0e8, alpha=1.2e-5)
        sec_truss = Section(id="truss_sec", A=0.01, I=0.0, d=0.0)
        truss = Element(id="T1", type="truss", node_i=self.node_i, node_j=self.node_j,
                       material=mat_truss, section=sec_truss)
        
        # Uniform 20°C temperature increase
        thermal = TemperatureL(element=truss, Tu=20.0, Tb=20.0)
        fef = thermal.FEF("pin-pin", 5.0)
        
        alpha = 1.2e-5
        E = 2.0e8
        A = 0.01
        T_uniform = 20.0
        F_T = alpha * T_uniform * E * A
        
        # Verify axial forces
        self.assertAlmostEqual(fef[0][0], -F_T, places=6)
        self.assertAlmostEqual(fef[3][0], F_T, places=6)
        # Verify no transverse/moment effects
        self.assertAlmostEqual(fef[1][0], 0.0, places=10)
        self.assertAlmostEqual(fef[2][0], 0.0, places=10)

    def test_thermal_gradient_temperature_frame(self):
        """TEST 2: Temperature gradient in frame element.
        
        Verifies combined axial and bending effects.
        Top = +10°C, Bottom = -10°C (delta_T = -20°C)
        Formula: M_T = (alpha * delta_T / d) * E * I
        """
        frame = Element(id="B1", type="frame", node_i=self.node_i, node_j=self.node_j,
                       material=self.mat_frame, section=self.sec_frame)
        
        # Gradient: top warmer
        thermal = TemperatureL(element=frame, Tu=10.0, Tb=-10.0)
        fef = thermal.FEF("fixed-fixed", 8.0)
        
        delta_T = -10.0 - 10.0  # Tb - Tu = -20
        T_uniform = 10.0 + (delta_T / 2.0)  # = 0
        
        alpha = 1.0e-5
        E = 3.0e10
        A = 0.18
        I = 0.0054
        d = 0.6
        
        F_T = alpha * T_uniform * E * A
        M_T = (alpha * delta_T / d) * E * I
        
        # Verify: [-F_T, 0, -M_T, F_T, 0, M_T]^T
        self.assertAlmostEqual(fef[0][0], -F_T, places=4)
        self.assertAlmostEqual(fef[2][0], -M_T, places=4)
        self.assertAlmostEqual(fef[3][0], F_T, places=4)
        self.assertAlmostEqual(fef[5][0], M_T, places=4)


class TestSupportSettlement(unittest.TestCase):
    """Unit tests for support settlement constraints."""
    
    def setUp(self):
        self.node_1 = Node(id=1, x=0.0, y=0.0)
        self.node_2 = Node(id=2, x=5.0, y=0.0)
        self.node_3 = Node(id=3, x=10.0, y=0.0)

    def test_settlement_horizontal_displacement_only(self):
        """TEST 1: Support with horizontal settlement only.
        
        Verifies: Support node constrained with UX settlement = -0.01 m
        Expected: Node has restrain_ux=True and settlement_ux=-0.01
        """
        support = Support(node=self.node_1,
                         restrain_ux=True,
                         restrain_uy=True,
                         restrain_rz=True,
                         settlement_ux=-0.01,
                         settlement_uy=0.0,
                         settlement_rz=0.0)
        
        self.assertTrue(support.restrain_ux)
        self.assertEqual(support.settlement_ux, -0.01)
        self.assertEqual(support.settlement_uy, 0.0)
        self.assertEqual(support.settlement_rz, 0.0)

    def test_settlement_vertical_and_rotation(self):
        """TEST 2: Support with vertical settlement and rotation.
        
        Verifies: Support node constrained with UY settlement = -0.005 m and RZ = -0.001 rad
        Expected: Node has restrain_uy=True, settlement_uy=-0.005, settlement_rz=-0.001
        """
        support = Support(node=self.node_2,
                         restrain_ux=False,
                         restrain_uy=True,
                         restrain_rz=True,
                         settlement_ux=0.0,
                         settlement_uy=-0.005,
                         settlement_rz=-0.001)
        
        self.assertFalse(support.restrain_ux)
        self.assertTrue(support.restrain_uy)
        self.assertTrue(support.restrain_rz)
        self.assertEqual(support.settlement_uy, -0.005)
        self.assertEqual(support.settlement_rz, -0.001)


if __name__ == '__main__':
    unittest.main()