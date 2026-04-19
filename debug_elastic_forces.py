#!/usr/bin/env python3
"""Debug script to trace K*u calculation for T1 element"""

import sys
import os
sys.path.insert(0, os.path.abspath('src'))

from parser import XMLParser
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver
from post_processor import PostProcessor
from element_physics import ElementPhysics
import math_utils

xml_path = "data/Assignment_4_Q2b.xml"
model = XMLParser(xml_path).parse()

# Get T1 element
t1 = model.elements['T1']
et1 = ElementPhysics(t1)

# Run full analysis to get displacements
opt = DOFOptimizer(model)
num_eq, semi_bw, _ = opt.optimize()
assembler = MatrixAssembler(model, num_eq, semi_bw)
K_banded, F_global = assembler.assemble("LC1")
solver = BandedSolver(K_banded, F_global, semi_bw)
D_active = solver.solve()
processor = PostProcessor(model, D_active, "LC1")

# Extract displacements for T1 nodes
u_node_i = [processor.displacements[t1.node_i.id][0], processor.displacements[t1.node_i.id][1]]
u_node_j = [processor.displacements[t1.node_j.id][0], processor.displacements[t1.node_j.id][1]]

print(f"T1 Element Displacements:")
print(f"  Node I (Global): UX={u_node_i[0]:.9e} m, UY={u_node_i[1]:.9e} m")
print(f"  Node J (Global): UX={u_node_j[0]:.9e} m, UY={u_node_j[1]:.9e} m")
print()

# Transform to local displacements
c = (t1.node_j.x - t1.node_i.x) / et1.L
s = (t1.node_j.y - t1.node_i.y) / et1.L

u_local_i_x = c * u_node_i[0] + s * u_node_i[1]
u_local_i_y = -s * u_node_i[0] + c * u_node_i[1]

u_local_j_x = c * u_node_j[0] + s * u_node_j[1]
u_local_j_y = -s * u_node_j[0] + c * u_node_j[1]

print(f"T1 Element Displacements (Local):")
print(f"  Node I: u_x={u_local_i_x:.9e} m, u_y={u_local_i_y:.9e} m")
print(f"  Node J: u_x={u_local_j_x:.9e} m, u_y={u_local_j_y:.9e} m")
print()

# Get local stiffness
k_local = et1.get_local_k()
print(f"Local Stiffness Matrix (4x4 truss):")
for i in range(4):
    print(f"  K[{i}] = [{k_local[i][0]:12.3e}, {k_local[i][1]:12.3e}, {k_local[i][2]:12.3e}, {k_local[i][3]:12.3e}]")
print()

# Calculate elastic force from K*u
u_local = [[u_local_i_x], [u_local_i_y], [u_local_j_x], [u_local_j_y]]
f_elastic = math_utils.matmul(k_local, u_local)

print(f"Elastic Forces (K*u):")
for i in range(4):
    print(f"  f_elastic[{i}] = {f_elastic[i][0]:12.6f} kN")
print()

# Get FEF
lc = model.load_cases['LC1']
fef_local = et1.get_local_fef(lc, model)

print(f"Fixed-End Forces (FEF):")
for i in range(4):
    print(f"  fef_local[{i}] = {fef_local[i][0]:12.6f} kN")
print()

# Total forces
f_total = math_utils.add(f_elastic, fef_local)

print(f"Total Member Forces (K*u + FEF):")
for i in range(4):
    print(f"  f_total[{i}] = {f_total[i][0]:12.6f} kN")
print()

print(f"Member forces from post-processor:")
forces = processor.member_forces['T1']
for i in range(4):
    print(f"  forces[{i}] = {forces[i][0]:12.6f} kN")
