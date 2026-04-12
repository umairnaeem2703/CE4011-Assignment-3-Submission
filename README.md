# Structural Analysis Solver (Python)

A matrix-based 2D structural analysis project for **truss**, **frame**, and **beam** models defined in XML.  
The solver parses input data, optimizes DOF numbering for bandwidth reduction, assembles a banded global system, solves it, and writes engineering result reports.

---

## Overview

The main workflow is implemented in [`run_analysis`](src/main.py) inside [src/main.py](src/main.py).  
At a high level, each run performs:

1. XML parsing  
2. DOF optimization / equation numbering  
3. Banded stiffness matrix assembly  
4. Linear solve  
5. Post-processing and report generation

---

## Key Features

- XML-driven model definition (nodes, elements, supports, load cases)
- Supports multiple structural problem types in one workspace
  - Truss-style models
  - Frame-style models
  - Beam-style models
- Per-load-case result output with:
  - Nodal displacements (`UX`, `UY`, `RZ`)
  - Member local end forces (Axial, Shear, Moment)
  - Support reactions (`Fx`, `Fy`, `Mz`)
- Batch-ready architecture via reusable analysis function: [`run_analysis`](src/main.py)

---

## Repository Layout

- Solver entry point: [src/main.py](src/main.py)
- Core solver package: [src/](src/)
- Input datasets: [data/](data/)
- Generated reports: [results/](results/)
- Tests: [tests/](tests/)

---

## How to Run

From workspace root:

```bash
python src/main.py
```

By default, [`main`](src/main.py) runs `./data/example1_case1_truss.xml` if it exists.

To analyze a different XML model programmatically, call [`run_analysis`](src/main.py) from [src/main.py](src/main.py).

---

## Inputs

Example input files are available in [data/](data/), including:

- [data/example1_case1_truss.xml](data/example1_case1_truss.xml)
- [data/example2_frame.xml](data/example2_frame.xml)
- [data/example3_frame_truss.xml](data/example3_frame_truss.xml)
- [data/Q3a.xml](data/Q3a.xml)
- [data/Q3b.xml](data/Q3b.xml)
- [data/Q3c.xml](data/Q3c.xml)
- [data/Q3d.xml](data/Q3d.xml)

---

## Outputs

Analysis reports are written to [results/](results/) using the naming pattern:

`{ModelName}_{LoadCase}_results.txt`

Examples:

- [results/Example1_Case1_Truss_LC1_results.txt](results/Example1_Case1_Truss_LC1_results.txt)
- [results/Q3c_L_Frame_Cantilever_LC1_results.txt](results/Q3c_L_Frame_Cantilever_LC1_results.txt)
- [results/Q3d_A_Frame_Pinned_LC1_results.txt](results/Q3d_A_Frame_Pinned_LC1_results.txt)
- [results/Q3e_Beam_All_Pins_LC1_results.txt](results/Q3e_Beam_All_Pins_LC1_results.txt)

Each report contains:

1. **Nodal Displacements**
2. **Member Local End Forces**
3. **Support Reactions**

---

## Assignment Utilities

These are useful for generating or validating input structures separate from the main solver flow in [src/main.py](src/main.py).

---

## Notes


- Keep XML inputs in [data/](data/) and use [results/](results/) for generated reports.
- If parsing or stability issues occur, check constraints/support conditions and connectivity in your XML model first.

---

## Quick Start Checklist

- Place or choose an XML file in [data/](data/).
- Run `python src/main.py`.
- Open generated report in [results/](results/).