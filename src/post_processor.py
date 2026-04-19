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
        """Maps the solved active equations back to all nodes (including 0.0 for restrained DOFs)."""
        for n_id, node in self.model.nodes.items():
            disp = []
            for dof in node.dofs:
                if dof >= 0:
                    disp.append(self.D_active[dof][0])
                else:
                    disp.append(0.0)
            self.displacements[n_id] = disp

    def _compute_forces_and_reactions(self):
        """Calculates local member forces and global support reactions."""
        # Initialize reaction sums at supports
        for n_id in self.model.supports.keys():
            self.reactions[n_id] = [0.0, 0.0, 0.0]

        for el_id, el in self.model.elements.items():
            phys = ElementPhysics(el)
            
            # 1. Retrieve condensed matrices (pass model for auto FEF detection)
            k_local = phys.get_local_k()
            fef_local = phys.get_local_fef(self.load_case, self.model)
            k_cond, fef_cond = phys.condense(k_local, fef_local)
            
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
                
            # 3. Calculate Local Forces: {f'} = [k_cond]{d'} + {FEF_cond}
            d_local = math_utils.matmul(T, d_global)
            f_local = math_utils.add(math_utils.matmul(k_cond, d_local), fef_cond)
            self.member_forces[el_id] = f_local
            
            # 4. Transform to global forces and sum at supports for reactions
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
                    
        # 5. Subtract applied nodal point loads from the calculated reactions
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
                
                f_i = [forces[0][0], 0.0, 0.0] if el.type == 'truss' else [forces[0][0], forces[1][0], forces[2][0]]
                f_j = [forces[2][0], 0.0, 0.0] if el.type == 'truss' else [forces[3][0], forces[4][0], forces[5][0]]
                
                f.write(f"{el_id:<8} {el.node_i.id:<4} {f_i[0]:<16.4f} {f_i[1]:<16.4f} {f_i[2]:<16.4f}\n")
                f.write(f"{'':<8} {el.node_j.id:<4} {f_j[0]:<16.4f} {f_j[1]:<16.4f} {f_j[2]:<16.4f}\n")
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