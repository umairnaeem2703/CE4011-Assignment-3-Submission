# src/main.py

import sys
import os
from parser import XMLParser
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver, UnstableStructureError
from post_processor import PostProcessor

def run_analysis(xml_filepath: str, output_dir: str = "./results"):
    print(f"--- Starting Analysis: {os.path.basename(xml_filepath)} ---")
    
    # 1. Parse Data
    try:
        parser = XMLParser(xml_filepath)
        model = parser.parse()
    except Exception as e:
        print(f"❌ Parse Error: {e}")
        return

    # 2. Optimize DOFs
    try:
        optimizer = DOFOptimizer(model)
        num_eq, semi_bw, full_bw = optimizer.optimize()
    except UnstableStructureError as e:
        print(f"❌ Structural Error:\n{e}")
        return

    print(f"Nodes: {len(model.nodes)} | Elements: {len(model.elements)}")
    print(f"Active Equations: {num_eq} | Bandwidth: {full_bw}")

    os.makedirs(output_dir, exist_ok=True)

    # Process each Load Case
    for lc_id in model.load_cases.keys():
        print(f"\nProcessing Load Case: {lc_id}...")
        
        # 3. Assemble Global System
        assembler = MatrixAssembler(model, num_eq, semi_bw)
        K_banded, F_global = assembler.assemble(lc_id)
        
        # 4. Solve System
        solver = BandedSolver(K_banded, F_global, semi_bw)
        try:
            D_active = solver.solve()
            print(f"✅ System solved successfully.")
        except UnstableStructureError as e:
            print(f"❌ {e}")
            continue
        except Exception as e:
            print(f"❌ Solver failed: {e}")
            continue

        # 5. Post-Process & Output
        processor = PostProcessor(model, D_active, lc_id)
        output_file = os.path.join(output_dir, f"{model.name}_{lc_id}_results.txt")
        processor.write_results(output_file)

if __name__ == "__main__":
    # Test on one of your provided files
    test_file = "./data/example3_frame_truss.xml"
    if os.path.exists(test_file):
        run_analysis(test_file)
    else:
        print("Please provide a valid path to an XML file in the script.")