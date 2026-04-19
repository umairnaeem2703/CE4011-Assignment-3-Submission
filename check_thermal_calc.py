# Quick thermal FEF calculation check

# For T1 element: uniform thermal load delta_T = 50°C
# - alpha (steel) = 1.2e-5 (1/°C)
# - E (steel) = 200000000 (Pa or units in input)
# - A (truss) = 0.000091106 (m²)
# - Tu = 50 (from uniform load parsing)
# - Tb = 50 (from uniform load parsing)

Tu = 50.0
Tb = 50.0
delta_T = Tb - Tu  # = 0
T_uniform = Tu + (delta_T / 2.0)  # = 50

alpha = 1.2e-5
E = 200000000
A = 0.000091106

F_T = alpha * T_uniform * E * A

print(f"Tu = {Tu}, Tb = {Tb}")
print(f"delta_T (Tb - Tu) = {delta_T}")
print(f"T_uniform = {T_uniform}")
print(f"F_T = {alpha} * {T_uniform} * {E} * {A}")
print(f"F_T = {F_T} kN")
print()
print("SAP2000 Reference: T1 Axial = -4.067 kN (at node 2)")
print(f"Ratio SAP/Computed: {4.067 / max(F_T, 0.01)}")
print()

# Check if maybe the gradient is being computed wrong
# The uniform load type parser sets Tu=dT, Tb=dT
# But maybe it should be set differently?

print("Alternative interpretation:")
print("If 'uniform' means delta_T is the change from reference (0):")
Tu_alt = 0.0
Tb_alt = 50.0
delta_T_alt = Tb_alt - Tu_alt
T_uniform_alt = Tu_alt + (delta_T_alt / 2.0)
F_T_alt = alpha * T_uniform_alt * E * A
print(f"Tu = {Tu_alt}, Tb = {Tb_alt}")
print(f"T_uniform = {T_uniform_alt}")
print(f"F_T = {F_T_alt} kN")
print()

print("Check: what T would give F = 4.067 kN?")
T_needed = 4.067 / (alpha * E * A)
print(f"Needed T = {T_needed} °C")
