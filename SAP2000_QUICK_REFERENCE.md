# SAP2000 Test Utility - Quick Reference

## One-Minute Overview

### The Three Essential Transformations

#### 1. Coordinate System Mapping (Automatic in Parser)
```
SAP2000 (X-Z plane)  →  Our Solver (X-Y plane)
U1, U3, R2 (SAP)     →  UX, UY, RZ (ours)
F1, F3, M2 (SAP)     →  Fx, Fy, Mz (ours)

Sign transforms:
  RZ = -R2  (Rotation sign flip)
  Mz = -M2  (Moment sign flip)
```

#### 2. Local-to-Global Force Transformation
```python
For element at angle θ:
  Fx = cos(θ)*Axial - sin(θ)*Shear
  Fy = sin(θ)*Axial + cos(θ)*Shear
  Mz = Moment (unchanged)
```

#### 3. Tolerance Application
```
Displacements: rel_tol = 0.01    (1% - different algorithms)
Forces:        rel_tol = 0.02    (2% - shear deformation gap)
Rotations:     rel_tol = 0.05    (5% - bending theory coupling)
Absolute:      abs_tol = 1e-6    (1 μm floor for tiny values)
```

## Usage Pattern

```python
from sap2000_parser import SAP2000Parser, MemberForceTransformer, assert_displacement_match

# Step 1: Parse SAP2000 file (coordinate mapping applied automatically)
parser = SAP2000Parser("sap2000_file.txt")
sap_disp, sap_react, sap_forces = parser.parse()

# Step 2: Run solver
processor = solve_model(model)

# Step 3: Compare displacements (direct, no transformation needed)
for node_id in sap_disp:
    match, error = assert_displacement_match(
        processor.displacements[node_id],
        sap_disp[node_id],
        rel_tol=0.01
    )

# Step 4: Transform member forces and compare
theta, cos_theta, sin_theta, L = MemberForceTransformer.calculate_element_angle(
    element.node_i, element.node_j
)
our_global = MemberForceTransformer.local_to_global_forces(
    axial, shear, moment, cos_theta, sin_theta
)
# Compare our_global to sap_forces[elem_id]['i' or 'j']
```

## Classes and Methods

### SAP2000Parser
```python
parser = SAP2000Parser(filepath)
displacements, reactions, element_forces = parser.parse()
# Returns:
#   displacements: {node_id: (ux, uy, rz)}
#   reactions:     {node_id: (fx, fy, mz)}
#   element_forces: {elem_id: {'i': (fx, fy, mz), 'j': (fx, fy, mz)}}
```

### CoordinateMapper
```python
ux, uy, rz = CoordinateMapper.map_sap2000_displacement(u1, u3, r2)
fx, fy, mz = CoordinateMapper.map_sap2000_force(f1, f3, m2)
```

### MemberForceTransformer
```python
# Calculate element orientation
theta, cos_theta, sin_theta, length = MemberForceTransformer.calculate_element_angle(
    node_i, node_j
)

# Transform local to global
fx, fy, mz = MemberForceTransformer.local_to_global_forces(
    axial, shear, moment, cos_theta, sin_theta
)

# Transform global to local (inverse)
axial, shear, moment = MemberForceTransformer.global_to_local_forces(
    fx, fy, mz, cos_theta, sin_theta
)
```

### Comparison Helpers
```python
match, error_msg = assert_displacement_match(computed, expected, rel_tol=0.01)
match, error_msg = assert_force_match(computed, expected, rel_tol=0.02)

# Both return (bool, str) - match status and error details
```

## Common Test Pattern

```python
class TestRegression(unittest.TestCase):
    def setUp(self):
        self.disp_rel_tol = 0.01
        self.force_rel_tol = 0.02
    
    def test_with_sap2000(self):
        # Parse and run
        parser = SAP2000Parser(sap_path)
        sap_disp, sap_react, sap_forces = parser.parse()
        processor = run_analysis(model)
        
        # Validate displacements
        for node_id in sap_disp:
            match, error = assert_displacement_match(
                processor.displacements[node_id],
                sap_disp[node_id],
                rel_tol=self.disp_rel_tol
            )
            self.assertTrue(match, f"Node {node_id}: {error}")
        
        # Validate reactions  
        for node_id in sap_react:
            match, error = assert_force_match(
                processor.reactions[node_id],
                sap_react[node_id],
                rel_tol=self.force_rel_tol
            )
            self.assertTrue(match, f"Node {node_id}: {error}")
        
        # Validate member forces (with transformation)
        for elem_id in sap_forces:
            if len(processor.member_forces.get(elem_id, [])) < 6:
                continue  # Skip truss elements
            
            theta, cos_t, sin_t, _ = MemberForceTransformer.calculate_element_angle(
                model.elements[elem_id].node_i,
                model.elements[elem_id].node_j
            )
            
            our_global = MemberForceTransformer.local_to_global_forces(
                processor.member_forces[elem_id][0][0],  # P_i
                processor.member_forces[elem_id][1][0],  # V_i
                processor.member_forces[elem_id][2][0],  # M_i
                cos_t, sin_t
            )
            
            match, error = assert_force_match(
                our_global,
                sap_forces[elem_id]['i'],
                rel_tol=self.force_rel_tol
            )
            self.assertTrue(match, f"Element {elem_id}: {error}")
```

## Key Implementation Details

**Why -R2 for RZ?**
- SAP2000's R2 uses right-hand rule looking from +Y (into the X-Z plane)
- Our RZ uses right-hand rule perpendicular to X-Y plane
- Different axis orientation = sign flip

**Why transform local to global?**
- SAP2000 exports global member forces (F1, F3, M2)
- Our solver computes local element forces (Axial, Shear, Moment)
- Before comparison, transform ours to match SAP2000's reference frame

**Why hybrid tolerance?**
- Relative tolerance fails on very small values (1% of 1e-6 is 1e-8)
- Absolute tolerance unnecessary for large values
- Using max(relative_error ≤ rel_tol OR absolute_error ≤ abs_tol) handles both

## Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `IndexError: list index out of range` | Accessing force component [3-5] on truss | Check `len(forces) < 6` and skip |
| `AssertionError: displacement mismatch` | Tolerance too tight or coordinate error | Increase `rel_tol` or check mapping |
| `Wrong sign on forces` | Missed sign convention | Use `CoordinateMapper` class |
| `Member forces don't match` | Forgot local→global transform | Apply `local_to_global_forces()` |

## Files

- **Parser**: `tests/sap2000_parser.py` (480 lines, fully documented)
- **Tests**: `tests/test_regression.py` (updated tests 1-6)
- **Guide**: `SAP2000_TESTING_GUIDE.md` (detailed documentation)

