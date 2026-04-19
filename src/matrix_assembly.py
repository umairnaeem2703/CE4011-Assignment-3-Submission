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

        return K_banded, F_global