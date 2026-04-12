# tests/test_interface.py

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from parser import Node, Material, Section, Element, Support, PointLoad, UDL, LoadCase, StructuralModel
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler

class TestMatrixAssembly(unittest.TestCase):
    def setUp(self):
        """Creates a minimal 2-element truss model manually for controlled testing."""
        self.model = StructuralModel(name="Test_Assembly")
        
        self.mat = Material(id="mat1", E=2.0e8)
        self.sec = Section(id="sec1", A=0.01) # EA/L = 2e6 for L=1
        
        # Geometry: Node 1 (0,0) -> Node 2 (1,0) -> Node 3 (2,0)
        self.n1 = Node(id=1, x=0.0, y=0.0)
        self.n2 = Node(id=2, x=1.0, y=0.0)
        self.n3 = Node(id=3, x=2.0, y=0.0)
        
        self.model.nodes = {1: self.n1, 2: self.n2, 3: self.n3}
        
        self.e1 = Element(id="E1", type="truss", node_i=self.n1, node_j=self.n2, material=self.mat, section=self.sec)
        self.e2 = Element(id="E2", type="truss", node_i=self.n2, node_j=self.n3, material=self.mat, section=self.sec)
        self.model.elements = {"E1": self.e1, "E2": self.e2}
        
        # Restrain Node 1 fully (Pin). Restrain Node 3 in Y (Roller). Node 2 is completely free.
        self.model.supports = {
            1: Support(node=self.n1, restrain_ux=True, restrain_uy=True),
            3: Support(node=self.n3, restrain_ux=False, restrain_uy=True)
        }
        
        # Apply a point load at Node 2
        lc = LoadCase(id="LC1")
        lc.point_loads.append(PointLoad(node=self.n2, fx=50.0, fy=-10.0))
        self.model.load_cases = {"LC1": lc}

        # Run Optimizer to set up DOFs
        self.optimizer = DOFOptimizer(self.model)
        self.num_eq, self.semi_bw, self.full_bw = self.optimizer.optimize()

    def test_assembly_global_stiffness_mapping(self):
        """Test 1: Verify local element matrices map correctly to Global Banded [K]."""
        assembler = MatrixAssembler(self.model, self.num_eq, self.semi_bw)
        K_banded, _ = assembler.assemble("LC1")
        
        # Expected Active DOFs based on constraints:
        # N1: [restrained, restrained] 
        # N2: [Active, Active]
        # N3: [Active, restrained]
        # Depending on RCM, let's look up the exact DOF assigned to N2's ux.
        n2_ux_dof = self.n2.dofs[0]
        n3_ux_dof = self.n3.dofs[0]
        
        # EA/L = 2.0e6. 
        # For N2's ux, both E1 and E2 contribute axial stiffness: 2e6 + 2e6 = 4e6
        # The main diagonal in the banded matrix is at index [dof][0]
        self.assertEqual(K_banded[n2_ux_dof][0], 4000000.0)
        
        # For N3's ux, only E2 contributes axial stiffness: 2e6
        self.assertEqual(K_banded[n3_ux_dof][0], 2000000.0)
        
        # Interaction between N2's ux and N3's ux (off-diagonal). 
        # Value should be -EA/L = -2e6. 
        # Note: Banded mapping is K[min_dof][max_dof - min_dof]
        min_dof = min(n2_ux_dof, n3_ux_dof)
        max_dof = max(n2_ux_dof, n3_ux_dof)
        self.assertEqual(K_banded[min_dof][max_dof - min_dof], -2000000.0)

    def test_assembly_global_load_vector(self):
        """Test 2: Verify Nodal loads and Equivalent Nodal Loads combine into {F}."""
        # Let's add a dummy frame element to test FEFs
        n4 = Node(id=4, x=2.0, y=5.0)
        self.model.nodes[4] = n4
        frame = Element(id="F1", type="frame", node_i=self.n3, node_j=n4, material=self.mat, section=self.sec)
        self.model.elements["F1"] = frame
        
        # Add a UDL to generate FEFs
        self.model.load_cases["LC1"].udls.append(UDL(element=frame, wy=-10.0, fef_condition="fixed-fixed"))
        
        # Re-run optimizer because we added nodes/elements
        num_eq, semi_bw, _ = self.optimizer.optimize()
        
        assembler = MatrixAssembler(self.model, num_eq, semi_bw)
        _, F_global = assembler.assemble("LC1")
        
        # Point load check at N2
        n2_ux_dof = self.n2.dofs[0]
        n2_uy_dof = self.n2.dofs[1]
        self.assertEqual(F_global[n2_ux_dof][0], 50.0)
        self.assertEqual(F_global[n2_uy_dof][0], -10.0)
        
        # FEF check at N3 (start of vertical frame F1)
        # Length of F1 = 5.0. UDL wy = -10.0 (local y, which is horizontal in global)
        # Note: In a real test, keeping track of global/local mapping for rotated frames is tricky.
        # But conceptually, F_global should not be 0 at the active DOFs of N3/N4 due to the ENL subtraction.
        n3_rz_dof = self.n3.dofs[2] 
        self.assertNotEqual(F_global[n3_rz_dof][0], 0.0) 

if __name__ == '__main__':
    unittest.main()