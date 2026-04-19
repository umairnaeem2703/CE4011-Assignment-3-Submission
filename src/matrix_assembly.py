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
        Assembles banded stiffness matrix [K] and load vector {F}.
        Handles element contributions, nodal loads, and settlement forces via matrix partitioning.
        """
        safe_bandwidth = max(self.semi_bandwidth + 1, self.num_active_dofs)
        K_banded = math_utils.zeros(self.num_active_dofs, safe_bandwidth)
        F_global = math_utils.zeros(self.num_active_dofs, 1)

        load_case = self.model.load_cases.get(load_case_id)
        if not load_case:
            raise ValueError(f"Load case '{load_case_id}' not found in the model.")

        # Element stiffness and fixed-end force contributions
        for element_id, element in self.model.elements.items():
            physics = ElementPhysics(element)
            
            k_local = physics.get_local_k()
            fef_local = physics.get_local_fef(load_case, self.model)
            
            k_condensed, fef_condensed = physics.condense(k_local, fef_local)
            k_global, fef_global = physics.transform_to_global(k_condensed, fef_condensed)
            
            if element.type == 'truss':
                element_dofs = element.node_i.dofs[0:2] + element.node_j.dofs[0:2]
            else:
                element_dofs = element.node_i.dofs + element.node_j.dofs
            
            for row_idx in range(len(element_dofs)):
                dof_row = element_dofs[row_idx]
                
                if dof_row >= 0:
                    F_global[dof_row][0] -= fef_global[row_idx][0]
                    
                    for col_idx in range(len(element_dofs)):
                        dof_col = element_dofs[col_idx]
                        
                        if dof_col >= dof_row:
                            band_col = dof_col - dof_row
                            K_banded[dof_row][band_col] += k_global[row_idx][col_idx]

        # Nodal load contributions
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

        # Settlement forces via matrix partitioning
        # Compute element-level unbalanced forces: {f} = [k]_e * {u_prescribed}
        # where u_prescribed contains settlements at restrained DOFs
        for element_id, element in self.model.elements.items():
            if element.type == 'truss':
                element_dofs = element.node_i.dofs[0:2] + element.node_j.dofs[0:2]
                num_dofs = 4
            else:
                element_dofs = element.node_i.dofs + element.node_j.dofs
                num_dofs = 6
            
            u_prescribed = math_utils.zeros(num_dofs, 1)
            has_settlement = False
            
            # Assemble prescribed displacements from supports
            node_i_support = self.model.supports.get(element.node_i.id)
            if node_i_support is not None:
                if node_i_support.restrain_ux:
                    u_prescribed[0][0] = node_i_support.settlement_ux
                    has_settlement = has_settlement or (node_i_support.settlement_ux != 0.0)
                if node_i_support.restrain_uy:
                    u_prescribed[1][0] = node_i_support.settlement_uy
                    has_settlement = has_settlement or (node_i_support.settlement_uy != 0.0)
            
            node_j_offset = 2 if element.type == 'truss' else 3
            node_j_support = self.model.supports.get(element.node_j.id)
            if node_j_support is not None:
                if node_j_support.restrain_ux:
                    u_prescribed[node_j_offset + 0][0] = node_j_support.settlement_ux
                    has_settlement = has_settlement or (node_j_support.settlement_ux != 0.0)
                if node_j_support.restrain_uy:
                    u_prescribed[node_j_offset + 1][0] = node_j_support.settlement_uy
                    has_settlement = has_settlement or (node_j_support.settlement_uy != 0.0)
            
            if not has_settlement:
                continue
            
            # Compute element unbalanced forces and subtract from F_global
            physics = ElementPhysics(element)
            k_local = physics.get_local_k()
            fef_local = physics.get_local_fef(load_case, self.model)
            k_condensed, _ = physics.condense(k_local, fef_local)
            k_global, _ = physics.transform_to_global(k_condensed, math_utils.zeros(num_dofs, 1))
            
            f_unbalanced_global = math_utils.matmul(k_global, u_prescribed)
            
            for dof_idx in range(len(element_dofs)):
                dof = element_dofs[dof_idx]
                if dof >= 0:
                    F_global[dof][0] -= f_unbalanced_global[dof_idx][0]

        return K_banded, F_global