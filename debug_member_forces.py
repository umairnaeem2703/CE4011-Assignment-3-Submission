#!/usr/bin/env python3
"""Debug script to trace member forces during post-processing"""

import sys
import os
sys.path.insert(0, os.path.abspath('src'))

from parser import XMLParser
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver
from post_processor import PostProcessor

xml_path = "data/Assignment_4_Q2b.xml"
model = XMLParser(xml_path).parse()

# Run the full analysis
opt = DOFOptimizer(model)
num_eq, semi_bw, _ = opt.optimize()

assembler = MatrixAssembler(model, num_eq, semi_bw)
K_banded, F_global = assembler.assemble("LC1")

solver = BandedSolver(K_banded, F_global, semi_bw)
D_active = solver.solve()

processor = PostProcessor(model, D_active, "LC1")

# Check node displacements
print("Node Displacements:")
for node_id in sorted(processor.displacements.keys()):
    disp = processor.displacements[node_id]
    print(f"  Node {node_id}: UX={disp[0]:.9e}, UY={disp[1]:.9e}, RZ={disp[2]:.9e}")

# Check member forces
if 'T1' in processor.member_forces:
    print("T1 Member Forces (from post-processor):")
    forces = processor.member_forces['T1']
    print(f"  Force at node I (local x): f[0] = {forces[0][0]:.6f} kN")
    print(f"  Force at node I (local y): f[1] = {forces[1][0]:.6f} kN")
    print(f"  Force at node J (local x): f[2] = {forces[2][0]:.6f} kN")
    print(f"  Force at node J (local y): f[3] = {forces[3][0]:.6f} kN")
else:
    print("T1 not in member_forces!")
    print(f"Available members: {list(processor.member_forces.keys())}")

# Also check T2 to see the pattern
if 'T2' in processor.member_forces:
    print("\nT2 Member Forces (from post-processor):")
    forces = processor.member_forces['T2']
    print(f"  Force at node I (local x): f[0] = {forces[0][0]:.6f} kN")
    print(f"  Force at node J (local x): f[2] = {forces[2][0]:.6f} kN")
    print("\nExpected from SAP2000:")
    print(f"  T2 node I: 5.676 kN")
    print(f"  T2 node J: -5.676 kN")
