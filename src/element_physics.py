# src/element_physics.py

import math
import math_utils
from parser import Element, LoadCase, MemberLoad, StructuralModel

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

    def _determine_fef_condition(self, model: StructuralModel) -> str:
        """Auto-detects the correct FEF condition based on element releases and global supports."""
        if self.element.type == 'truss':
            return "pin-pin"

        rel_i = self.element.release_start
        rel_j = self.element.release_end

        if model is not None:
            # Check if node_i is a boundary pin with no other frame connections
            sup_i = model.supports.get(self.element.node_i.id)
            if sup_i and not sup_i.restrain_rz:
                frames_at_i = sum(1 for e in model.elements.values() 
                                  if e.type == 'frame' and 
                                  (e.node_i.id == self.element.node_i.id or e.node_j.id == self.element.node_i.id))
                if frames_at_i == 1:
                    rel_i = True

            # Check if node_j is a boundary pin with no other frame connections
            sup_j = model.supports.get(self.element.node_j.id)
            if sup_j and not sup_j.restrain_rz:
                frames_at_j = sum(1 for e in model.elements.values() 
                                  if e.type == 'frame' and 
                                  (e.node_i.id == self.element.node_j.id or e.node_j.id == self.element.node_j.id))
                if frames_at_j == 1:
                    rel_j = True

        if rel_i and rel_j: return "pin-pin"
        if rel_i: return "pin-fixed"
        if rel_j: return "fixed-pin"
        return "fixed-fixed"

    def get_local_fef(self, load_case: LoadCase, model: StructuralModel = None) -> list:
        """Calculates combined Fixed End Forces (FEF) delegated to the load subclasses."""
        if self.element.type == 'truss':
            fef_total = math_utils.zeros(4, 1)
            fef_cond = "pin-pin"  # Trusses always use pin-pin for FEF
        else:
            fef_total = math_utils.zeros(6, 1)
            fef_cond = self._determine_fef_condition(model)

        for load in load_case.loads:
            if isinstance(load, MemberLoad) and load.element.id == self.element.id:
                # For temperature loads, always use "fixed-fixed" since thermal stresses are built-in
                # regardless of boundary conditions
                load_fef_cond = "fixed-fixed" if hasattr(load, 'Tu') else fef_cond
                
                # Load instances calculate their own FEF based on the appropriate condition
                fef_member = load.FEF(load_fef_cond, self.L)
                
                # For trusses, extract only the first 4 components (axial only)
                if self.element.type == 'truss':
                    fef_member = [[fef_member[i][0]] for i in range(4)]
                
                fef_total = math_utils.add(fef_total, fef_member)

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