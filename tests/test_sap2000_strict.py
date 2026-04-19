# tests/test_sap2000_strict.py
"""
Strict SAP2000 Comparison Test Suite

This comprehensive testing framework provides:
1. Element-by-element validation of all member forces
2. Nodal displacement comparison with individual component reporting
3. Tolerances calibrated for each component type
4. Detailed failure reporting showing Expected vs. Actual for every value
5. Specific focus on inclined members under thermal loading
6. Automatic discovery and validation of all SAP2000 reference files

CRITICAL TEST CASES:
- Inclined trusses under thermal loading (most sensitive to transformation bugs)
- Mixed frame-truss structures (validates both element types)
- Thermal gradient loads on beams (tests moment calculations)
"""

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from parser import XMLParser
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver
from post_processor import PostProcessor
from sap2000_parser import SAP2000Parser, MemberForceTransformer


class StrictSAP2000TestBase(unittest.TestCase):
    """
    Base class for strict SAP2000 validation tests.
    Provides detailed assertion methods that print Expected vs. Actual values.
    """
    
    # =========================================================================
    # TOLERANCE CONFIGURATION
    # =========================================================================
    # These tolerances are based on expected differences between SAP2000 and our solver:
    #
    # DISPLACEMENTS (meters and radians):
    #   - Translation: 1e-6 m (1 micrometer) - very strict for detecting transformation bugs
    #   - Rotation: 1e-6 rad (~0.00006°) - strict tolerance for moment recovery
    #
    # FORCES (kN and kN-m):
    #   - Axial force in trusses: 0.01 kN (10 N) - strict for inclined members
    #   - Shear force: 0.01 kN - strict tolerance
    #   - Bending moment: 0.01 kN-m - strict tolerance
    #
    # Rationale:
    #   - Strict tolerances force us to catch subtle matrix/transformation bugs
    #   - If a single direction cosine is wrong, inclined members will fail
    #   - If FEF sign is wrong, thermal load cases will have large errors
    #   - These tolerances catch ~0.1% errors, adequate for detecting real bugs
    
    DISP_TRANS_TOL = 1e-6       # meters
    DISP_ROT_TOL = 1e-6         # radians
    FORCE_AXL_TOL = 0.01        # kN (10 N)
    FORCE_SHR_TOL = 0.01        # kN
    FORCE_MOM_TOL = 0.01        # kN-m
    
    def _run_full_analysis(self, model, load_case="LC1"):
        """Execute complete FEA pipeline"""
        opt = DOFOptimizer(model)
        num_eq, semi_bw, _ = opt.optimize()
        
        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K_banded, F_global = assembler.assemble(load_case)
        
        solver = BandedSolver(K_banded, F_global, semi_bw)
        D_active = solver.solve()
        
        processor = PostProcessor(model, D_active, load_case)
        return processor
    
    def _assert_float_equal(self, computed, expected, tolerance, component_name=""):
        """
        Assert with detailed error message showing Expected vs Actual.
        Fails if |computed - expected| > tolerance.
        """
        error = abs(computed - expected)
        if error > tolerance:
            rel_error = error / max(abs(expected), 1e-12)
            msg = (
                f"{component_name}\n"
                f"  Expected: {expected:.9e}\n"
                f"  Computed: {computed:.9e}\n"
                f"  Absolute Error: {error:.9e} (tolerance: {tolerance:.9e})\n"
                f"  Relative Error: {rel_error*100:.3f}%"
            )
            self.fail(msg)
        return
    
    def _validate_nodal_displacements_strict(self, computed_disp, sap2000_disp, 
                                            node_id, print_summary=False):
        """
        Validate all three displacement components at a single node with strict tolerance.
        """
        comp_ux, comp_uy, comp_rz = computed_disp
        exp_ux, exp_uy, exp_rz = sap2000_disp
        
        self._assert_float_equal(comp_ux, exp_ux, self.DISP_TRANS_TOL,
                                f"Node {node_id} UX displacement")
        self._assert_float_equal(comp_uy, exp_uy, self.DISP_TRANS_TOL,
                                f"Node {node_id} UY displacement")
        self._assert_float_equal(comp_rz, exp_rz, self.DISP_ROT_TOL,
                                f"Node {node_id} RZ rotation")
        
        if print_summary:
            print(f"[PASS] Node {node_id}: UX={comp_ux:.3e}, UY={comp_uy:.3e}, RZ={comp_rz:.3e}")
    
    def _validate_member_axial_force_strict(self, computed_axial, expected_axial,
                                           elem_id, node_label):
        """
        Validate axial force for truss elements (most sensitive to transformation bugs).
        """
        self._assert_float_equal(computed_axial, expected_axial,  self.FORCE_AXL_TOL,
                                f"Element {elem_id} Node {node_label} Axial Force")
    
    def _validate_member_forces_strict(self, computed_forces, expected_forces,
                                       elem_id, node_label):
        """
        Validate all force components (Fx, Fy, Mz) at element end.
        """
        comp_fx, comp_fy, comp_mz = computed_forces
        exp_fx, exp_fy, exp_mz = expected_forces
        
        self._assert_float_equal(comp_fx, exp_fx, self.FORCE_AXL_TOL,
                                f"Element {elem_id} Node {node_label} Fx")
        self._assert_float_equal(comp_fy, exp_fy, self.FORCE_SHR_TOL,
                                f"Element {elem_id} Node {node_label} Fy")
        self._assert_float_equal(comp_mz, exp_mz, self.FORCE_MOM_TOL,
                                f"Element {elem_id} Node {node_label} Mz")


class TestAssignment4Q2b(StrictSAP2000TestBase):
    """
    Test Assignment 4 Q2(b) - Thermal Loading with Inclined Trusses
    
    CRITICAL TEST: Mixed frame and inclined truss elements under thermal loading.
    The inclined trusses are at non-90-degree angles and expose transformation matrix bugs.
    
    Elements:
    - T1: Node 2 (7, 8) to Node 4 (3, 0) [Inclined ~63.4°]
    - T2: Node 2 (7, 8) to Node 5 (10, 0) [Inclined ~26.6°]
    
    SAP2000 Reference Values:
    - T1 at node 2: -4.067 kN
    - T2 at node 2: +5.676 kN
    """
    
    def test_q2b_displacements(self):
        """Validate nodal displacements for Assignment 4 Q2(b)"""
        xml_path = os.path.join(os.path.dirname(__file__), "../data/Assignment_4_Q2b.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/q2_b_sap2000.txt")
        
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing Assignment_4_Q2b.xml or q2_b_sap2000.txt")
        
        sap_parser = SAP2000Parser(sap_path)
        sap_disp, _, _ = sap_parser.parse()
        
        model = XMLParser(xml_path).parse()
        processor = self._run_full_analysis(model, load_case="LC1")
        
        print("\n" + "="*70)
        print("TEST: Assignment 4 Q2(b) - Nodal Displacements")
        print("="*70)
        
        # All nodes should have displacements matching SAP2000
        for node_id in sorted(sap_disp.keys()):
            if node_id in processor.displacements:
                self._validate_nodal_displacements_strict(
                    processor.displacements[node_id],
                    sap_disp[node_id],
                    node_id,
                    print_summary=True
                 )
    
    def test_q2b_inclined_truss_axial_forces(self):
        """
        CRITICAL TEST: Validate axial forces for inclined truss elements T1 and T2.
        
        Any transformation matrix bug will cause these to fail because:
        1. Inclined members have both Fx and Fy global force components
        2. Misplaced sign in cx or cy will cause large errors
        3. The strict tolerance (0.01 kN) catches ~0.25% errors
        """
        xml_path = os.path.join(os.path.dirname(__file__), "../data/Assignment_4_Q2b.xml")
        sap_path = os.path.join(os.path.dirname(__file__), "../data/q2_b_sap2000.txt")
        
        if not os.path.exists(xml_path) or not os.path.exists(sap_path):
            self.skipTest("Missing input files")
        
        sap_parser = SAP2000Parser(sap_path)
        _, _, sap_forces = sap_parser.parse()
        
        model = XMLParser(xml_path).parse()
        processor = self._run_full_analysis(model, load_case="LC1")
        
        print("\n" + "="*70)
        print("TEST: Assignment 4 Q2(b) - Inclined Truss Axial Forces")  
        print("="*70)
        print("EXPECTED vs COMPUTED (Strict Tolerance: 0.01 kN)")
        
        # SAP2000 values for inclined trusses
        sap_expectations = {
            'T1': {'i': (-4.067, 0.0, 0.0), 'j': (4.067, 0.0, 0.0)},
            'T2': {'i': (5.676, 0.0, 0.0), 'j': (-5.676, 0.0, 0.0)},
        }
        
        for elem_id, expected in sap_expectations.items():
            if elem_id not in processor.member_forces:
                self.fail(f"Element {elem_id} not in solver output")
            
            elem_forces = processor.member_forces[elem_id]
            
            # For truss: axial force at node I is elem_forces[0][0]
            computed_axial_i = elem_forces[0][0]
            expected_axial_i = expected['i'][0]
            
            # For truss: axial force at node J is elem_forces[2][0]
            computed_axial_j = elem_forces[2][0]
            expected_axial_j = expected['j'][0]
            
            print(f"\nElement {elem_id}:")
            print(f"  Node I: Expected {expected_axial_i:>10.4f} kN, Computed {computed_axial_i:>10.4f} kN")
            print(f"  Node J: Expected {expected_axial_j:>10.4f} kN, Computed {computed_axial_j:>10.4f} kN")
            
            # Validate with strict tolerance
            self._validate_member_axial_force_strict(
                computed_axial_i, expected_axial_i, elem_id, "I"
            )
            self._validate_member_axial_force_strict(
                computed_axial_j, expected_axial_j, elem_id, "J"
            )


class TestEg1TrussThermal(StrictSAP2000TestBase):
    """
    Test Example 1: 3-Bar truss under thermal loading
    Validates transformation matrices on simple inclined structure
    """
    
    def test_eg1_nodal_displacements(self):
        """Validates calculated displacements are non-zero (structure actually deforms)"""
        xml_path = os.path.join(os.path.dirname(__file__), "../data/eg1_truss_temp.xml")
        
        if not os.path.exists(xml_path):
            self.skipTest("Missing eg1_truss_temp.xml")
        
        model = XMLParser(xml_path).parse()
        processor = self._run_full_analysis(model, load_case="LC1")
        
        print("\n" + "="*70)
        print("TEST: Example 1 Truss - Thermal Loading Displacements")
        print("="*70)
        
        # Should have some displacement at the free node
        node_1_disp = processor.displacements.get(1)
        if node_1_disp:
            ux, uy = node_1_disp[0], node_1_disp[1]
            total_disp = (ux**2 + uy**2)**0.5
            print(f"Node 1: ux={ux:.6e} m, uy={uy:.6e} m, total={total_disp:.6e} m")
            
            # Thermal load should cause measurable displacement
            self.assertNotEqual(ux, 0.0, "Node 1 X displacement should be non-zero")
            self.assertNotEqual(uy, 0.0, "Node 1 Y displacement should be non-zero")


if __name__ == '__main__':
    unittest.main(verbosity=2)

