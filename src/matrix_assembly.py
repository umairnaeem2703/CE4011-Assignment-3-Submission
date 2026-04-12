# src/matrix_assembly.py

import math_utils
from parser import StructuralModel
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
        # Initialize the restricted banded matrix and load vector
        # Allocate bandwidth generously to handle all cases safely
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
            
            # 1a. Get local matrices
            k_local = physics.get_local_k()
            fef_local = physics.get_local_fef(load_case)
            
            # 1b. Apply static condensation for internal hinges
            k_condensed, fef_condensed = physics.condense(k_local, fef_local)
            
            # 1c. Transform to global coordinates
            k_global, fef_global = physics.transform_to_global(k_condensed, fef_condensed)
            
            # 1d. Map to Global Banded [K] and {F}
            # Gather DOFs for this element based on element type
            if element.type == 'truss':
                # Truss uses only 2 DOFs per node (ux, uy)
                element_dofs = element.node_i.dofs[0:2] + element.node_j.dofs[0:2]
            else:
                # Frame uses all 3 DOFs per node (ux, uy, rz)
                element_dofs = element.node_i.dofs + element.node_j.dofs
            
            for row_idx in range(len(element_dofs)):
                dof_row = element_dofs[row_idx]
                
                # Skip restrained DOFs (-1)
                if dof_row >= 0:
                    # Subtract FEF from global force vector ({F} = {P} - {FEF})
                    F_global[dof_row][0] -= fef_global[row_idx][0]
                    
                    # Add to banded stiffness matrix
                    for col_idx in range(len(element_dofs)):
                        dof_col = element_dofs[col_idx]
                        
                        # We only store the upper triangular part for the banded solver
                        if dof_col >= dof_row:
                            # The column index in the banded array is (dof_col - dof_row)
                            band_col = dof_col - dof_row
                            
                            K_banded[dof_row][band_col] += k_global[row_idx][col_idx]

        # ==========================================
        # 2. ASSEMBLE NODAL POINT LOADS
        # ==========================================
        for point_load in load_case.point_loads:
            node = point_load.node
            
            # ux
            if node.dofs[0] >= 0:
                F_global[node.dofs[0]][0] += point_load.fx
            # uy
            if node.dofs[1] >= 0:
                F_global[node.dofs[1]][0] += point_load.fy
            # rz (only applicable to frame nodes)
            if len(node.dofs) > 2 and node.dofs[2] >= 0:
                F_global[node.dofs[2]][0] += point_load.mz

        return K_banded, F_global