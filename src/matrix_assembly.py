# src/matrix_assembly.py

import math_utils
from parser import StructuralModel, NodalLoad
from element_physics import ElementPhysics

class MatrixAssembler:
    def __init__(self, model: StructuralModel, num_active_dofs: int, semi_bandwidth: int):
        self.model = model
        self.num_active_dofs = num_active_dofs
        self.semi_bandwidth = semi_bandwidth

    def assemble(self, load_case_id: str) -> tuple[list, list]:
        """
        Assembles the banded global stiffness matrix [K] and load vector {F}.
        Implements matrix partitioning for support settlements:
          - K is partitioned into K_ff (free × free), K_fr (free × restrained)
          - For prescribed settlements U_r at restrained DOFs:
            K_ff * U_f = F - K_fr * U_r
          - Settlement forces are computed element-by-element without forming K_fr
        
        Returns: K_banded (2D list), F_global (2D list as column vector).
        """
        safe_bandwidth = max(self.semi_bandwidth + 1, self.num_active_dofs)
        K_banded = math_utils.zeros(self.num_active_dofs, safe_bandwidth)
        F_global = math_utils.zeros(self.num_active_dofs, 1)

        load_case = self.model.load_cases.get(load_case_id)
        if not load_case:
            raise ValueError(f"Load case '{load_case_id}' not found in the model.")

        # ==========================================
        # 1. ASSEMBLE ELEMENT CONTRIBUTIONS
        # ==========================================
        for element_id, element in self.model.elements.items():
            physics = ElementPhysics(element)
            
            k_local = physics.get_local_k()
            # Pass model so physics can auto-determine FEF condition
            fef_local = physics.get_local_fef(load_case, self.model)
            
            k_condensed, fef_condensed = physics.condense(k_local, fef_local)
            k_global, fef_global = physics.transform_to_global(k_condensed, fef_condensed)
            
            if element.type == 'truss':
                element_dofs = element.node_i.dofs[0:2] + element.node_j.dofs[0:2]
            else:
                element_dofs = element.node_i.dofs + element.node_j.dofs
            
            # Assemble K_banded and F_global contributions
            for row_idx in range(len(element_dofs)):
                dof_row = element_dofs[row_idx]
                
                if dof_row >= 0:
                    F_global[dof_row][0] -= fef_global[row_idx][0]
                    
                    for col_idx in range(len(element_dofs)):
                        dof_col = element_dofs[col_idx]
                        
                        if dof_col >= dof_row:
                            band_col = dof_col - dof_row
                            K_banded[dof_row][band_col] += k_global[row_idx][col_idx]

        # ==========================================
        # 2. ASSEMBLE NODAL LOADS (Directly)
        # ==========================================
        for load in load_case.loads:
            if isinstance(load, NodalLoad):
                node = load.node
                nodal_forces = load.NodalLoads()
                
                if node.dofs[0] >= 0:
                    F_global[node.dofs[0]][0] += nodal_forces[0]
                if node.dofs[1] >= 0:
                    F_global[node.dofs[1]][0] += nodal_forces[1]
                if len(node.dofs) > 2 and node.dofs[2] >= 0:
                    F_global[node.dofs[2]][0] += nodal_forces[2]

        # ==========================================
        # 3. APPLY SETTLEMENT FORCES (Matrix Partitioning)
        # ==========================================
        # For each element, check if any of its nodes have non-zero settlements.
        # Compute element-level unbalanced forces: {f_unbalanced}_e = [k]_e * {u_r}_e
        # where u_r contains settlements at restrained DOFs (0 at free DOFs).
        # Subtract the active DOF components from F_global.
        # This avoids forming the global K_fr matrix while achieving:
        #   F_modified = F_original - K_fr * U_r
        
        for element_id, element in self.model.elements.items():
            # Build the known displacement vector for this element
            # (settlements at restrained DOFs, zeros elsewhere)
            if element.type == 'truss':
                element_dofs = element.node_i.dofs[0:2] + element.node_j.dofs[0:2]
                num_dofs = 4
            else:
                element_dofs = element.node_i.dofs + element.node_j.dofs
                num_dofs = 6
            
            # Assemble known displacements for this element
            # (only non-zero where settlements exist at restrained DOFs)
            u_prescribed = math_utils.zeros(num_dofs, 1)
            has_settlement = False
            
            # Node i displacements
            node_i_support = self.model.supports.get(element.node_i.id)
            if node_i_support is not None:
                if node_i_support.restrain_ux:
                    u_prescribed[0][0] = node_i_support.settlement_ux
                    has_settlement = has_settlement or (node_i_support.settlement_ux != 0.0)
                if node_i_support.restrain_uy:
                    u_prescribed[1][0] = node_i_support.settlement_uy
                    has_settlement = has_settlement or (node_i_support.settlement_uy != 0.0)
                if num_dofs == 6 and node_i_support.restrain_rz:
                    # Note: rotation settlements can be added here if needed
                    pass
            
            # Node j displacements
            node_j_offset = 2 if element.type == 'truss' else 3
            node_j_support = self.model.supports.get(element.node_j.id)
            if node_j_support is not None:
                if node_j_support.restrain_ux:
                    u_prescribed[node_j_offset + 0][0] = node_j_support.settlement_ux
                    has_settlement = has_settlement or (node_j_support.settlement_ux != 0.0)
                if node_j_support.restrain_uy:
                    u_prescribed[node_j_offset + 1][0] = node_j_support.settlement_uy
                    has_settlement = has_settlement or (node_j_support.settlement_uy != 0.0)
                if num_dofs == 6 and node_j_support.restrain_rz:
                    # Note: rotation settlements can be added here if needed
                    pass
            
            # Skip if no settlements for this element
            if not has_settlement:
                continue
            
            # Compute element unbalanced forces: {f_unbalanced} = [k]_e * {u_prescribed}
            physics = ElementPhysics(element)
            k_local = physics.get_local_k()
            fef_local = physics.get_local_fef(load_case, self.model)
            k_condensed, _ = physics.condense(k_local, fef_local)
            k_global, _ = physics.transform_to_global(k_condensed, math_utils.zeros(num_dofs, 1))
            
            # {f_unbalanced}_e = [k]_e * {u_prescribed}_e
            f_unbalanced_global = math_utils.matmul(k_global, u_prescribed)
            
            # Subtract settlement forces from active DOF load vector
            # Only DOFs that are active (dof >= 0) are affected
            for dof_idx in range(len(element_dofs)):
                dof = element_dofs[dof_idx]
                if dof >= 0:
                    # Subtract unbalanced force from this active DOF
                    F_global[dof][0] -= f_unbalanced_global[dof_idx][0]

        return K_banded, F_global