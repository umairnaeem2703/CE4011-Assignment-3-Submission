# src/post_processor.py

import os
import math_utils
from element_physics import ElementPhysics
from parser import StructuralModel, NodalLoad

class PostProcessor:
    def __init__(self, model: StructuralModel, D_active: list, load_case_id: str):
        self.model = model
        self.D_active = D_active
        self.load_case = model.load_cases[load_case_id]
        
        self.displacements = {}  # node_id -> [ux, uy, rz]
        self.member_forces = {}  # el_id -> list of local forces
        self.reactions = {}      # node_id -> [Fx, Fy, Mz]
        
        self._build_full_displacements()
        self._compute_forces_and_reactions()

    def _build_full_displacements(self):
        """
        Maps the solved active DOF displacements back to all nodes,
        including prescribed settlement values at restrained DOFs.
        
        For each node:
          - Active DOFs: use the solved value from D_active
          - Restrained DOFs with settlement: use the prescribed settlement value
          - Restrained DOFs without settlement: use 0.0
        """
        for n_id, node in self.model.nodes.items():
            disp = []
            for dof_idx, dof in enumerate(node.dofs):
                if dof >= 0:
                    # Active DOF: use solved value
                    disp.append(self.D_active[dof][0])
                else:
                    # Restrained DOF: check for prescribed settlement
                    support = self.model.supports.get(n_id)
                    if support is not None:
                        if dof_idx == 0 and support.restrain_ux:
                            # X-direction restrained with possible settlement
                            disp.append(support.settlement_ux)
                        elif dof_idx == 1 and support.restrain_uy:
                            # Y-direction restrained with possible settlement
                            disp.append(support.settlement_uy)
                        elif dof_idx == 2 and support.restrain_rz:
                            # Rotation restrained (settlement_rz not yet implemented)
                            disp.append(0.0)
                        else:
                            disp.append(0.0)
                    else:
                        # No support at this node, restrained DOF = 0.0
                        disp.append(0.0)
            self.displacements[n_id] = disp

    def _compute_forces_and_reactions(self):
        """
        Calculates local member forces and global support reactions using the correct pipeline:
        1. Extract global displacements u_global (4x1 for truss, 6x1 for frame)
        2. Transform to local: u_local = T @ u_global
        3. Compute member forces: f_local = (k_local @ u_local) + fef_local
        4. Extract axial/shear/moment from the full local force vector
        """
        # Initialize reaction sums at supports
        for n_id in self.model.supports.keys():
            self.reactions[n_id] = [0.0, 0.0, 0.0]

        for el_id, el in self.model.elements.items():
            phys = ElementPhysics(el)
            
            # 1. Get full local stiffness and FEF (not condensed)
            k_local = phys.get_local_k()
            fef_local = phys.get_local_fef(self.load_case, self.model)
            
            # 2. Build rotation matrix [T] and extract global displacements
            c, s = phys.cos_x, phys.sin_x
            if el.type == 'truss':
                T = [
                    [ c,  s,  0,  0],
                    [-s,  c,  0,  0],
                    [ 0,  0,  c,  s],
                    [ 0,  0, -s,  c]
                ]
                d_global = [
                    [self.displacements[el.node_i.id][0]], [self.displacements[el.node_i.id][1]],
                    [self.displacements[el.node_j.id][0]], [self.displacements[el.node_j.id][1]]
                ]
            else:
                T = [
                    [ c,  s,  0,  0,  0,  0], [-s,  c,  0,  0,  0,  0], [ 0,  0,  1,  0,  0,  0],
                    [ 0,  0,  0,  c,  s,  0], [ 0,  0,  0, -s,  c,  0], [ 0,  0,  0,  0,  0,  1]
                ]
                d_global = [
                    [self.displacements[el.node_i.id][0]], [self.displacements[el.node_i.id][1]], [self.displacements[el.node_i.id][2]],
                    [self.displacements[el.node_j.id][0]], [self.displacements[el.node_j.id][1]], [self.displacements[el.node_j.id][2]]
                ]
            
            # 3. Transform displacements to local coordinates: u_local = T @ u_global
            u_local = math_utils.matmul(T, d_global)
            
            # 4. Compute local member forces: f_local = (k_local @ u_local) + fef_local
            f_local = math_utils.add(math_utils.matmul(k_local, u_local), fef_local)
            self.member_forces[el_id] = f_local
            
            # 5. Transform local forces back to global for reaction calculation
            f_global = math_utils.matmul(math_utils.transpose(T), f_local)
            
            if el.node_i.id in self.reactions:
                self.reactions[el.node_i.id][0] += f_global[0][0]
                self.reactions[el.node_i.id][1] += f_global[1][0]
                if el.type == 'frame': self.reactions[el.node_i.id][2] += f_global[2][0]
                    
            if el.node_j.id in self.reactions:
                offset = 2 if el.type == 'truss' else 3
                self.reactions[el.node_j.id][0] += f_global[offset + 0][0]
                self.reactions[el.node_j.id][1] += f_global[offset + 1][0]
                if el.type == 'frame': self.reactions[el.node_j.id][2] += f_global[offset + 2][0]
                    
        # 6. Subtract applied nodal point loads from the calculated reactions
        for load in self.load_case.loads:
            if isinstance(load, NodalLoad):
                n_id = load.node.id
                if n_id in self.reactions:
                    nodal_forces = load.NodalLoads()
                    self.reactions[n_id][0] -= nodal_forces[0]
                    self.reactions[n_id][1] -= nodal_forces[1]
                    self.reactions[n_id][2] -= nodal_forces[2]

    def write_results(self, filepath: str):
        """Formats and writes the engineering results to a text file."""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"STRUCTURAL ANALYSIS REPORT\n")
            f.write(f"Model: {self.model.name}\n")
            f.write(f"Load Case: {self.load_case.name} ({self.load_case.id})\n")
            f.write("=" * 60 + "\n\n")
            
            # --- NODAL DISPLACEMENTS ---
            f.write("1. NODAL DISPLACEMENTS\n")
            f.write("-" * 60 + "\n")
            f.write(f"{'Node':<8} {'UX (m)':<15} {'UY (m)':<15} {'RZ (rad)':<15}\n")
            for n_id in sorted(self.displacements.keys()):
                disp = self.displacements[n_id]
                f.write(f"{n_id:<8} {disp[0]:<15.6e} {disp[1]:<15.6e} {disp[2]:<15.6e}\n")
            f.write("\n")
            
            # --- MEMBER FORCES ---
            f.write("2. MEMBER LOCAL END FORCES\n")
            f.write("-" * 80 + "\n")
            f.write("Element  Node  Axial (kN)       Shear (kN)       Moment (kN-m)\n")
            f.write("-" * 80 + "\n")
            for el_id in sorted(self.member_forces.keys()):
                el = self.model.elements[el_id]
                forces = self.member_forces[el_id]
                
                if el.type == 'truss':
                    # Truss forces: [Fx_i, Fy_i, Fx_j, Fy_j]
                    f_i_axial = forces[0][0]
                    f_i_shear = 0.0
                    f_i_moment = 0.0
                    f_j_axial = forces[2][0]  # Fx at node J
                    f_j_shear = 0.0
                    f_j_moment = 0.0
                else:
                    # Frame forces: [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]
                    f_i_axial = forces[0][0]
                    f_i_shear = forces[1][0]
                    f_i_moment = forces[2][0]
                    f_j_axial = forces[3][0]  # Fx at node J (tension is positive)
                    f_j_shear = forces[4][0]
                    f_j_moment = forces[5][0]
                
                f.write(f"{el_id:<8} {el.node_i.id:<4} {f_i_axial:<16.4f} {f_i_shear:<16.4f} {f_i_moment:<16.4f}\n")
                f.write(f"{'':<8} {el.node_j.id:<4} {f_j_axial:<16.4f} {f_j_shear:<16.4f} {f_j_moment:<16.4f}\n")
                f.write("-" * 80 + "\n")
            f.write("\n")
            
            # --- SUPPORT REACTIONS ---
            f.write("3. SUPPORT REACTIONS\n")
            f.write("-" * 60 + "\n")
            f.write(f"{'Node':<8} {'Fx (kN)':<15} {'Fy (kN)':<15} {'Mz (kN-m)':<15}\n")
            for n_id in sorted(self.reactions.keys()):
                r = self.reactions[n_id]
                sup = self.model.supports[n_id]
                
                # Only report values if the DOF is actually restrained
                rx = f"{r[0]:.4f}" if sup.restrain_ux else "Free"
                ry = f"{r[1]:.4f}" if sup.restrain_uy else "Free"
                rz = f"{r[2]:.4f}" if sup.restrain_rz else "Free"
                
                f.write(f"{n_id:<8} {rx:<15} {ry:<15} {rz:<15}\n")

        print(f"✅ Results written to {filepath}")