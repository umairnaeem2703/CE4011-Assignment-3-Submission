#!/usr/bin/env python3
"""Debug script to trace actual FEF values being computed"""

import sys
sys.path.insert(0, 'src')

from parser import XMLParser
from element_physics import ElementPhysics

xml_path = "data/Assignment_4_Q2b.xml"
model = XMLParser(xml_path).parse()

# Get T1 element
t1 = model.elements['T1']
print(f"T1 Element: {t1.id}")
print(f"  Node I: {t1.node_i.id} ({t1.node_i.x}, {t1.node_i.y})")
print(f"  Node J: {t1.node_j.id} ({t1.node_j.x}, {t1.node_j.y})")
print(f"  Material: E={t1.material.E}, alpha={t1.material.alpha}")
print(f"  Section: A={t1.section.A}")
print()

# Get the load case
lc = model.load_cases['LC1']
print(f"Load case: {lc.id} - {lc.name}")
print(f"Loads: {len(lc.loads)}")
for load in lc.loads:
    if hasattr(load, 'element') and load.element.id == 'T1':
        print(f"  T1 Load: {type(load).__name__}")
        if hasattr(load, 'Tu'):
            print(f"    Tu={load.Tu}, Tb={load.Tb}")
print()

# Create element physics and calculate FEF
elem_phys = ElementPhysics(t1)
print(f"Element Physics:")
print(f"  Length L = {elem_phys.L:.4f}")
print(f"  Type: {elem_phys.element.type}")
print()

# Calculate FEF for the thermal load
for load in lc.loads:
    if hasattr(load, 'element') and load.element.id == 'T1':
        if hasattr(load, 'Tu'):
            fef_local = elem_phys.calculate_thermal_fef(load.Tu, load.Tb)
            
            # Manual calculation
            alpha = t1.material.alpha
            E = t1.material.E
            A = t1.section.A
            T_uniform = load.Tu + (load.Tb - load.Tu) / 2.0
            delta_T = load.Tb - load.Tu
            
            F_T = alpha * T_uniform * E * A
            
            print(f"Thermal Load Details:")
            print(f"  Tu={load.Tu}, Tb={load.Tb}")
            print(f"  delta_T = {delta_T}")  
            print(f"  T_uniform = {T_uniform}")
            print(f"  F_T = alpha * T_uniform * E * A")
            print(f"       = {alpha} * {T_uniform} * {E} * {A}")
            print(f"       = {F_T} kN")
            print()
            print(f"Calculated FEF_local:")
            for i, f in enumerate(fef_local):
                print(f"  [{i}] = {f[0]:.6f}")
            print()
            print(f"Expected FEF_local (from formula):")
            print(f"  [0] = {-F_T:.6f}  (compression at node I)")
            print(f"  [1] = 0.000000")
            print(f"  [2] = {F_T:.6f}  (compression at node J)")
            print(f"  [3] = 0.000000")
