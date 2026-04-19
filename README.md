# Structural Analysis Solver (Python)

A matrix-based 2D structural analysis solver for **truss**, **frame**, and **beam** models defined in XML.  
The program parses XML input files, optimizes node numbering for bandwidth reduction, assembles a banded system matrix, solves the linear system, and generates engineering reports with displacement and force results.

---

## Quick Start

### Prerequisites
- Python 3.10 or later
- Virtual environment (recommended)

### Setup & Execution

1. **Activate your Python virtual environment:**
   ```bash
   .venv\Scripts\activate
   ```

2. **Run the analysis:**
   ```bash
   python src/main.py
   ```

   This will analyze the default test file: `./data/Assignment_4_Q2b.xml`

3. **Check the results:**
   ```bash
   View the generated report in: ./results/
   ```

---

## Analysis Workflow

The solver executes the following steps for each load case:

1. **Parse XML** – Load structural model definition (nodes, elements, constraints, loads)
2. **Optimize DOF Numbering** – Reorder equations to minimize bandwidth
3. **Assemble System Matrix** – Build banded global stiffness matrix and load vector
4. **Solve Linear System** – Compute nodal displacements using banded solver
5. **Post-Process** – Calculate member forces and support reactions
6. **Generate Report** – Write results to output file

---

## Program Features

- **XML-Based Input** – Define models with nodes, elements, supports, and load cases
- **Multiple Problem Types** – Supports trusses, frames, and beams in one workspace
- **Comprehensive Output** – Per load case reports including:
  - Nodal displacements (UX, UY, RZ)
  - Member local end forces (Axial, Shear, Moment)
  - Support reactions (Fx, Fy, Mz)
- **Reusable API** – Call [`run_analysis()`](src/main.py) function for batch processing

---

## Project Structure

```
src/
├── main.py                    # Entry point & run_analysis() function
├── parser.py                  # XML model parsing
├── dof_optimizer.py          # DOF numbering optimization
├── matrix_assembly.py        # Stiffness matrix assembly
├── banded_solver.py          # Linear system solver
├── post_processor.py         # Results calculation & reporting
├── element_physics.py        # Element stiffness calculations
├── structural_validator.py   # Model validation
├── math_utils.py             # Linear algebra utilities
└── matrix_assembly.py        # Banded matrix assembly

data/                          # Input XML files
├── Assignment_4_Q2a.xml
├── Assignment_4_Q2b.xml
└── SCHEMA.xml               # XML schema documentation

results/                       # Generated output reports
tests/                        # Unit and integration tests
```

---

## Usage Examples

### Run Default Analysis
```bash
python src/main.py
```

### Analyze Specific XML File (Programmatically)
```python
from src.main import run_analysis

run_analysis("./data/Assignment_4_Q2a.xml", output_dir="./results")
```

### Run Tests
```bash
python -m pytest tests/
```

---

## Input Files

XML model files are located in [data/](data/). Each file defines:
- **Nodes** – Coordinates and boundary conditions (supports)
- **Elements** – Connectivity and material properties
- **Load Cases** – Nodal loads and distributed loads
- **Material** – Element stiffness properties

Example input files:
- `Assignment_4_Q2a.xml` – Truss analysis
- `Assignment_4_Q2b.xml` – Frame analysis

---

## Output Reports

Analysis results are written to [results/](results/) in the format:

```
{ModelName}_{LoadCase}_results.txt
```

Each report contains:
1. Model summary (nodes, elements, equations)
2. Nodal displacements table
3. Member local end forces table
4. Support reactions table

Example output:
- `Assignment_4_Q2a_LC1_results.txt`
- `Assignment_4_Q2b_LC1_results.txt`

---

## Core Functions

**[src/main.py](src/main.py) – `run_analysis(xml_filepath, output_dir)`**
- Executes a complete structural analysis
- Parameters:
  - `xml_filepath` (str) – Path to input XML file
  - `output_dir` (str) – Directory for output reports (default: `./results`)
- Returns: None (writes results to file)

---

## Error Handling

The program includes validation for:
- **Parse Errors** – Invalid XML or missing required fields
- **Unstable Structures** – Models with insufficient boundary conditions
- **Solver Failures** – Singular or ill-conditioned matrices

Error messages are printed to console with debugging context.

---

## Notes

- All output reports use SI units (meters, Newtons, Pascal)
- Bandwidth optimization reduces computation time for large models
- The solver uses direct banded Gaussian elimination
- Supports load case superposition for combined load analysis

---

## Notes


- Keep XML inputs in [data/](data/) and use [results/](results/) for generated reports.
- If parsing or stability issues occur, check constraints/support conditions and connectivity in your XML model first.

---

## Quick Start Checklist

- Place or choose an XML file in [data/](data/).
- Run `python src/main.py`.
- Open generated report in [results/](results/).