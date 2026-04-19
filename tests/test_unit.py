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

        w = -10.0
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
        
        P = -20.0
        a = 2.0
        L = 5.0
        b = L - a
        
        pt_load = PointLoad(element=frame, position=a, fy=P)
        fef_local = pt_load.FEF("fixed-fixed", L)
        
        expected_vy_i = (P * b**2 * (3*a + b)) / (L**3)
        expected_mz_i = (P * a * b**2) / (L**2)
        
        self.assertAlmostEqual(fef_local[1][0], expected_vy_i, places=3)
        self.assertAlmostEqual(fef_local[2][0], expected_mz_i, places=3)

    def test_fef_thermal_gradient(self):
        """Test 5: Compute FEFs for a thermal gradient load (Case 2)."""
        frame = Element(id="F1", type="frame", node_i=self.node_i, node_j=self.node_j, 
                        material=self.mat, section=self.sec)
        
        delta_T = 20.0 # Top warmer than bottom
        thermal_load = TemperatureL(element=frame, Tu=delta_T/2.0, Tb=-delta_T/2.0)
        fef_local = thermal_load.FEF("fixed-fixed", 5.0)

        # M_zi = -(alpha * delta_T / d) * EI
        expected_mz_i = -(1.2e-5 * 20.0 / 0.3) * (2.0e8 * 0.0001)
        
        self.assertAlmostEqual(fef_local[2][0], expected_mz_i, places=3)
        self.assertAlmostEqual(fef_local[5][0], -expected_mz_i, places=3)

if __name__ == '__main__':
    unittest.main()