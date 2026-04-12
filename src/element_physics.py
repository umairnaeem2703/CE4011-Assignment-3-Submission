# src/element_physics.py

import math
import math_utils
from parser import Element, LoadCase

class ElementPhysics:
    def __init__(self, element: Element):
        self.element = element
        self.L, self.cos_x, self.sin_x = self._calculate_geometry()

    def _calculate_geometry(self) -> tuple[float, float, float]:
        """Calculates length and direction cosines."""
        dx = self.element.node_j.x - self.element.node_i.x
        dy = self.element.node_j.y - self.element.node_i.y
        L = math.hypot(dx, dy)
        if L == 0:
            raise ValueError(f"Element {self.element.id} has zero length.")
        return L, dx / L, dy / L

    def get_local_k(self) -> list:
        """Returns the fully-fixed local stiffness matrix [k]."""
        E = self.element.material.E
        A = self.element.section.A
        L = self.L

        if self.element.type == 'truss':
            k = E * A / L
            return [
                [ k,  0, -k,  0],
                [ 0,  0,  0,  0],
                [-k,  0,  k,  0],
                [ 0,  0,  0,  0]
            ]
        else:
            I = self.element.section.I
            k_axial = E * A / L
            k_v1 = 12 * E * I / (L**3)
            k_v2 = 6 * E * I / (L**2)
            k_r1 = 4 * E * I / L
            k_r2 = 2 * E * I / L

            return [
                [ k_axial,       0,       0,-k_axial,       0,       0],
                [       0,    k_v1,    k_v2,       0,   -k_v1,    k_v2],
                [       0,    k_v2,    k_r1,       0,   -k_v2,    k_r2],
                [-k_axial,       0,       0, k_axial,       0,       0],
                [       0,   -k_v1,   -k_v2,       0,    k_v1,   -k_v2],
                [       0,    k_v2,    k_r2,       0,   -k_v2,    k_r1]
            ]

    def _compute_udl_fef(self, wy: float, wx: float, fef_condition: str) -> list:
        """Derives local Fixed End Forces (Equivalent Nodal Loads) for a UDL."""
        fef = math_utils.zeros(6, 1)
        L = self.L
        
        # Axial FEF (linear split)
        fef[0][0] = wx * L / 2.0
        fef[3][0] = wx * L / 2.0
        
        # Transverse FEF
        if fef_condition == "fixed-fixed":
            fef[1][0] = wy * L / 2.0
            fef[2][0] = wy * (L**2) / 12.0
            fef[4][0] = wy * L / 2.0
            fef[5][0] = -wy * (L**2) / 12.0
            
        elif fef_condition == "pin-fixed":
            fef[1][0] = (3.0 / 8.0) * wy * L
            fef[2][0] = 0.0
            fef[4][0] = (5.0 / 8.0) * wy * L
            fef[5][0] = wy * (L**2) / 8.0
            
        elif fef_condition == "fixed-pin":
            fef[1][0] = (5.0 / 8.0) * wy * L
            fef[2][0] = wy * (L**2) / 8.0
            fef[4][0] = (3.0 / 8.0) * wy * L
            fef[5][0] = 0.0
            
        elif fef_condition == "pin-pin":
            fef[1][0] = wy * L / 2.0
            fef[2][0] = 0.0
            fef[4][0] = wy * L / 2.0
            fef[5][0] = 0.0
            
        return fef

    def _compute_point_load_fef(self, fy: float, fx: float, a: float, fef_condition: str) -> list:
        """Derives local Fixed End Forces (Equivalent Nodal Loads) for a Member Point Load."""
        fef = math_utils.zeros(6, 1)
        L = self.L
        b = L - a
        
        # Axial FEF
        fef[0][0] = fx * b / L
        fef[3][0] = fx * a / L
        
        # Transverse FEF
        P = fy
        if fef_condition == "fixed-fixed":
            fef[1][0] = P * (b**2) * (3*a + b) / (L**3)
            fef[2][0] = P * a * (b**2) / (L**2)
            fef[4][0] = P * (a**2) * (3*b + a) / (L**3)
            fef[5][0] = -P * (a**2) * b / (L**2)
            
        elif fef_condition == "pin-fixed":
            Mz_j = P * a * b * (L + a) / (2 * L**2)
            fef[1][0] = (P * b - Mz_j) / L
            fef[2][0] = 0.0
            fef[4][0] = (P * a + Mz_j) / L
            fef[5][0] = Mz_j
            
        elif fef_condition == "fixed-pin":
            Mz_i = P * a * b * (L + b) / (2 * L**2)
            fef[1][0] = (P * b + Mz_i) / L
            fef[2][0] = Mz_i
            fef[4][0] = (P * a - Mz_i) / L
            fef[5][0] = 0.0
            
        elif fef_condition == "pin-pin":
            fef[1][0] = P * b / L
            fef[2][0] = 0.0
            fef[4][0] = P * a / L
            fef[5][0] = 0.0
            
        return fef

    def get_local_fef_point_load(self, P: float, a: float) -> list:
        """Direct wrapper strictly to support test_unit.py Test 5."""
        return self._compute_point_load_fef(fy=P, fx=0.0, a=a, fef_condition="fixed-fixed")

    def get_local_fef(self, load_case: LoadCase) -> list:
        """Calculates combined Fixed End Forces (FEF) from all member loads."""
        if self.element.type == 'truss':
            return math_utils.zeros(4, 1)

        fef_total = math_utils.zeros(6, 1)

        # Apply UDLs (Both v1 generic UDLs and v2 member_udl)
        for udl in load_case.udls:
            if udl.element.id == self.element.id:
                f_cond = getattr(udl, 'fef_condition', 'fixed-fixed')
                fef_udl = self._compute_udl_fef(udl.wy, udl.wx, f_cond)
                fef_total = math_utils.add(fef_total, fef_udl)

        # Apply point loads positioned along the member
        for mpl in load_case.member_point_loads:
            if mpl.element.id == self.element.id:
                f_cond = getattr(mpl, 'fef_condition', 'fixed-fixed')
                fef_mpl = self._compute_point_load_fef(mpl.fy, mpl.fx, mpl.position, f_cond)
                fef_total = math_utils.add(fef_total, fef_mpl)

        return fef_total

    def condense(self, k_local: list, fef_local: list) -> tuple[list, list]:
        """Performs generalized static condensation for hinged frames."""
        if self.element.type == 'truss':
            return k_local, fef_local 

        condensed_dofs = []
        if self.element.release_start:
            condensed_dofs.append(2) # Mz at node i
        if self.element.release_end:
            condensed_dofs.append(5) # Mz at node j

        if not condensed_dofs:
            return k_local, fef_local

        retained_dofs = [i for i in range(6) if i not in condensed_dofs]

        k_rr = [[k_local[r][c] for c in retained_dofs] for r in retained_dofs]
        k_cc = [[k_local[c_row][c_col] for c_col in condensed_dofs] for c_row in condensed_dofs]
        k_rc = [[k_local[r][c_col] for c_col in condensed_dofs] for r in retained_dofs]
        k_cr = [[k_local[c_row][c] for c in retained_dofs] for c_row in condensed_dofs]

        fef_r = [[fef_local[r][0]] for r in retained_dofs]
        fef_c = [[fef_local[c_row][0]] for c_row in condensed_dofs]

        k_cc_inv = math_utils.invert_matrix(k_cc)
        term1 = math_utils.matmul(math_utils.matmul(k_rc, k_cc_inv), k_cr)
        k_mod = math_utils.subtract(k_rr, term1)

        term2 = math_utils.matmul(math_utils.matmul(k_rc, k_cc_inv), fef_c)
        fef_mod = math_utils.subtract(fef_r, term2)

        k_final = math_utils.zeros(6, 6)
        fef_final = math_utils.zeros(6, 1)

        for i, r in enumerate(retained_dofs):
            fef_final[r][0] = fef_mod[i][0]
            for j, c in enumerate(retained_dofs):
                k_final[r][c] = k_mod[i][j]

        return k_final, fef_final

    def transform_to_global(self, k_local: list, fef_local: list) -> tuple[list, list]:
        """Rotates local matrices into the global coordinate system."""
        c, s = self.cos_x, self.sin_x
        
        if self.element.type == 'truss':
            T = [
                [ c,  s,  0,  0],
                [-s,  c,  0,  0],
                [ 0,  0,  c,  s],
                [ 0,  0, -s,  c]
            ]
        else:
            T = [
                [ c,  s,  0,  0,  0,  0],
                [-s,  c,  0,  0,  0,  0],
                [ 0,  0,  1,  0,  0,  0],
                [ 0,  0,  0,  c,  s,  0],
                [ 0,  0,  0, -s,  c,  0],
                [ 0,  0,  0,  0,  0,  1]
            ]

        T_trans = math_utils.transpose(T)
        
        k_global = math_utils.matmul(math_utils.matmul(T_trans, k_local), T)
        fef_global = math_utils.matmul(T_trans, fef_local)

        return k_global, fef_global

    def recover_local_forces(self, u_local: list, fef_local: list) -> list:
        """Recovers the local member-end forces: {F_local} = [k_local]{u_local} + {fef_local}."""
        k_local = self.get_local_k()
        ku = math_utils.matmul(k_local, u_local)
        return math_utils.add(ku, fef_local)