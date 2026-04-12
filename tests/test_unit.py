# tests/test_unit.py

import sys
import os
import unittest

# Ensure Python can find our source files
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from parser import Node, Material, Section, Element, LoadCase, UDL
from element_physics import ElementPhysics
import math_utils

class TestElementPhysics(unittest.TestCase):
    def setUp(self):
        # Setup dummy properties for testing
        self.mat = Material(id="steel", E=2.0e8)
        self.sec = Section(id="frame_sec", A=0.01, I=0.0001)
        # Default horizontal element (L=5) for standard tests
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

        # Verify specific matrix indices against theory (0-indexed)
        self.assertAlmostEqual(k_local[1][1], expected_k_v1, places=3)  # Shear stiffness
        self.assertAlmostEqual(k_local[2][2], expected_k_r1, places=3)  # Rotational stiffness

    def test_fef_uniform_distributed_load(self):
        """Test 2: Compute FEFs for a member with a UDL."""
        frame = Element(id="F1", type="frame", node_i=self.node_i, node_j=self.node_j, 
                        material=self.mat, section=self.sec)
        
        load_case = LoadCase(id="LC1")
        # Apply -10 kN/m transverse UDL
        load_case.udls.append(UDL(element=frame, wy=-10.0, fef_condition="fixed-fixed")) 
        
        physics = ElementPhysics(frame)
        fef_local = physics.get_local_fef(load_case)

        w = -10.0
        L = 5.0
        
        # Expected Vy = wL/2 (for fixed-fixed)
        expected_vy = (w * L) / 2.0
        # Expected Mz at start = wL^2 / 12 (for fixed-fixed)
        expected_mz_i = (w * L**2) / 12.0
        
        self.assertAlmostEqual(fef_local[1][0], expected_vy, places=3)
        self.assertAlmostEqual(fef_local[2][0], expected_mz_i, places=3)

    def test_transformation_matrix(self):
        """Test 3: Evaluate coordinate transformation for an inclined element."""
        # Create an inclined element: 3-4-5 triangle
        # length = 5, cos(theta) = 0.6, sin(theta) = 0.8
        node_inc_i = Node(id=1, x=0.0, y=0.0)
        node_inc_j = Node(id=2, x=3.0, y=4.0)
        
        frame = Element(id="F2", type="frame", node_i=node_inc_i, node_j=node_inc_j, 
                        material=self.mat, section=self.sec)
        physics = ElementPhysics(frame)
        
        # Verify direction cosines generated during init
        self.assertAlmostEqual(physics.cos_x, 0.6)
        self.assertAlmostEqual(physics.sin_x, 0.8)
        
        # Test 4b: Verify transform_to_global applies T matrix correctly
        # Let local FEF be purely axial: Fx_local = 10.0
        fef_local = [[10.0], [0.0], [0.0], [0.0], [0.0], [0.0]]
        k_local_dummy = math_utils.zeros(6, 6)
        
        _, fef_global = physics.transform_to_global(k_local_dummy, fef_local)
        
        # Global forces should be: Fx = 10*cos, Fy = 10*sin
        self.assertAlmostEqual(fef_global[0][0], 6.0)
        self.assertAlmostEqual(fef_global[1][0], 8.0)
        self.assertAlmostEqual(fef_global[2][0], 0.0)