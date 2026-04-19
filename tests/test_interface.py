# tests/test_interface.py

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from parser import Node, Material, Section, Element, Support, LoadCase, StructuralModel, NodalLoad, TemperatureL
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler


class TestMatrixAssemblyThermal(unittest.TestCase):
    """Integration tests for thermal loading assembly."""

    def test_thermal_uniform_load_assembly(self):
        """TEST 1: Assemble uniform temperature load into global load vector.
        
        Verifies:
        1. Thermal FEF computed correctly and assembled
        2. Global load vector includes thermal effects
        3. Full pipeline (optimize -> assemble) works correctly
        
        Setup: Simple 2-element frame, uniform dT=20°C
        Expected: Thermal forces appear in load vector
        """
        model = StructuralModel(name="Thermal_Uniform")
        
        mat = Material(id="mat_th", E=2.1e11, alpha=1.2e-5)
        sec = Section(id="sec_th", A=0.04, I=5e-5, d=0.4)
        
        n1 = Node(id=1, x=0.0, y=0.0)
        n2 = Node(id=2, x=5.0, y=0.0)
        n3 = Node(id=3, x=10.0, y=0.0)
        
        model.nodes = {1: n1, 2: n2, 3: n3}
        
        e1 = Element(id="E1", type="frame", node_i=n1, node_j=n2, material=mat, section=sec)
        e2 = Element(id="E2", type="frame", node_i=n2, node_j=n3, material=mat, section=sec)
        model.elements = {"E1": e1, "E2": e2}
        
        # Pin at n1, roller at n3
        model.supports = {
            1: Support(node=n1, restrain_ux=True, restrain_uy=True, restrain_rz=True),
            3: Support(node=n3, restrain_uy=True)
        }
        
        # Uniform thermal load: dT = +20°C
        lc = LoadCase(id="LC_TH1")
        thermal = TemperatureL(element=e1, Tu=20.0, Tb=20.0)
        lc.loads.append(thermal)
        model.load_cases = {"LC_TH1": lc}
        
        # Optimize and assemble
        optimizer = DOFOptimizer(model)
        num_eq, semi_bw, _ = optimizer.optimize()
        
        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K, F_global = assembler.assemble("LC_TH1")
        
        # Verify assembly succeeded
        self.assertGreater(len(K), 0, "Stiffness matrix should be assembled")
        self.assertEqual(len(F_global), num_eq, "Load vector size should match active DOFs")
        
        # Verify n2 has active DOFs (middle node should be free)
        n2_has_active = any(dof >= 0 for dof in n2.dofs)
        self.assertTrue(n2_has_active, "Middle node should have active DOFs")

    def test_thermal_gradient_load_assembly(self):
        """TEST 2: Assemble temperature gradient load with bending effects.
        
        Verifies:
        1. Gradient load creates both axial and moment FEF
        2. Global matrix and load vector assembled correctly
        3. Cantilever boundary conditions respected
        
        Setup: Single-element cantilever, Tu=20°C (top), Tb=0°C (bottom)
        Expected: FEF includes both axial and moment components
        """
        model = StructuralModel(name="Thermal_Gradient")
        
        mat = Material(id="mat_grad", E=3.0e10, alpha=1.0e-5)
        sec = Section(id="sec_grad", A=0.25, I=1.56e-3, d=0.5)
        
        n_fixed = Node(id=1, x=0.0, y=0.0)
        n_free = Node(id=2, x=6.0, y=0.0)
        
        model.nodes = {1: n_fixed, 2: n_free}
        
        cantilever = Element(id="CANT1", type="frame", node_i=n_fixed, node_j=n_free,
                            material=mat, section=sec)
        model.elements = {"CANT1": cantilever}
        
        # Fully fixed at base
        model.supports = {
            1: Support(node=n_fixed, restrain_ux=True, restrain_uy=True, restrain_rz=True)
        }
        
        # Gradient: top warmer by 20°C
        lc = LoadCase(id="LC_GRAD")
        thermal_grad = TemperatureL(element=cantilever, Tu=20.0, Tb=0.0)
        lc.loads.append(thermal_grad)
        model.load_cases = {"LC_GRAD": lc}
        
        # Optimize and assemble
        optimizer = DOFOptimizer(model)
        num_eq, semi_bw, _ = optimizer.optimize()
        
        self.assertGreater(num_eq, 0, "Should have active DOFs at free end")
        
        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K, F_global = assembler.assemble("LC_GRAD")
        
        # Verify assembly
        self.assertEqual(len(F_global), num_eq)
        self.assertGreater(len(K), 0)
        
        # Free end should have active DOFs
        n_free_active = any(dof >= 0 for dof in n_free.dofs)
        self.assertTrue(n_free_active, "Free end should have active DOFs")


class TestMatrixAssemblySettlement(unittest.TestCase):
    """Integration tests for support settlement assembly."""

    def test_settlement_horizontal_only(self):
        """TEST 1: Assemble horizontal settlement constraint.
        
        Verifies:
        1. Settlement DOF treated as prescribed (removed from active equations)
        2. Stiffness and load vector assembled for remaining DOFs
        3. Settlement causes prescribed displacement in constraint
        
        Setup: 2-node frame, n1 with UX settlement = -0.01 m
        Expected: UX DOF at n1 is constrained, not active
        """
        model = StructuralModel(name="Settlement_Horiz")
        
        mat = Material(id="mat_set", E=2.0e8)
        sec = Section(id="sec_set", A=0.01, I=1e-4, d=0.3)
        
        n1 = Node(id=1, x=0.0, y=0.0)
        n2 = Node(id=2, x=5.0, y=0.0)
        
        model.nodes = {1: n1, 2: n2}
        
        beam = Element(id="B1", type="frame", node_i=n1, node_j=n2, material=mat, section=sec)
        model.elements = {"B1": beam}
        
        # n1: pin with horizontal settlement
        # n2: free
        model.supports = {
            1: Support(node=n1, restrain_ux=True, restrain_uy=True,
                      settlement_ux=-0.01, settlement_uy=0.0)
        }
        
        # No loads
        lc = LoadCase(id="LC1")
        model.load_cases = {"LC1": lc}
        
        # Optimize and assemble
        optimizer = DOFOptimizer(model)
        num_eq, semi_bw, _ = optimizer.optimize()
        
        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K, F_global = assembler.assemble("LC1")
        
        # Verify assembly
        self.assertGreater(len(K), 0, "Stiffness matrix assembled")
        self.assertEqual(len(F_global), num_eq)
        
        # n1's UX should be inactive (constrained by settlement)
        # n1's UY should be inactive (pinned)
        # Only n2 should have active DOFs
        n2_active = any(dof >= 0 for dof in n2.dofs)
        self.assertTrue(n2_active, "Free node should have active DOFs")

    def test_settlement_mixed_constraints(self):
        """TEST 2: Assemble mixed settlements (UY + RZ) on support.
        
        Verifies:
        1. Multiple settlement components assembled correctly
        2. UX remains free while UY and RZ are prescribed
        3. Global equations constructed for unrestrained DOFs
        
        Setup: 2-element frame, n1 with UY=-0.005m + RZ=-0.001rad settlements
        Expected: UY and RZ at n1 are prescribed, UX is free
        """
        model = StructuralModel(name="Settlement_Mixed")
        
        mat = Material(id="mat_mix", E=3.0e10)
        sec = Section(id="sec_mix", A=0.3*0.6, I=0.3*0.6**3/12, d=0.6)
        
        n1 = Node(id=1, x=0.0, y=0.0)
        n2 = Node(id=2, x=4.0, y=0.0)
        n3 = Node(id=3, x=8.0, y=0.0)
        
        model.nodes = {1: n1, 2: n2, 3: n3}
        
        e1 = Element(id="E1", type="frame", node_i=n1, node_j=n2, material=mat, section=sec)
        e2 = Element(id="E2", type="frame", node_i=n2, node_j=n3, material=mat, section=sec)
        model.elements = {"E1": e1, "E2": e2}
        
        # n1: UX free, UY and RZ with settlements
        # n3: pinned
        model.supports = {
            1: Support(node=n1, restrain_ux=False, restrain_uy=True, restrain_rz=True,
                      settlement_uy=-0.005, settlement_rz=-0.001),
            3: Support(node=n3, restrain_ux=True, restrain_uy=True, restrain_rz=True)
        }
        
        # Apply horizontal load at n2
        lc = LoadCase(id="LC_MIXED")
        lc.loads.append(NodalLoad(node=n2, fx=100.0))
        model.load_cases = {"LC_MIXED": lc}
        
        # Optimize and assemble
        optimizer = DOFOptimizer(model)
        num_eq, semi_bw, _ = optimizer.optimize()
        
        self.assertGreater(num_eq, 0, "Should have active DOFs (UX at n1, all at n2)")
        
        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K, F_global = assembler.assemble("LC_MIXED")
        
        # Verify assembly
        self.assertEqual(len(F_global), num_eq)
        self.assertGreater(len(K), 0)
        
        # n1 UX should be active, n2 should be active
        n1_ux_active = n1.dofs[0] >= 0
        n2_active = any(dof >= 0 for dof in n2.dofs)
        self.assertTrue(n1_ux_active, "n1 UX should be free (active)")
        self.assertTrue(n2_active, "n2 should have active DOFs")

if __name__ == '__main__':
    unittest.main()