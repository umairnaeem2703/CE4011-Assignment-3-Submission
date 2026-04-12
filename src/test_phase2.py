# src/test_phase2.py

import os
from parser import XMLParser
from dof_optimizer import DOFOptimizer

def run_tests():
    # Adjust paths relative to the src/ directory
    test_files = [
        "./data/example1_case1_truss.xml",
        "./data/example1_case2_truss.xml",
        "./data/example2_frame.xml",
        "./data/example3_frame_truss.xml"
    ]

    print(f"{'Model Name':<35} | {'Active DOFs':<12} | {'Semi-BW (m)':<12} | {'Full-BW (b)':<12}")
    print("-" * 80)

    for filepath in test_files:
        if not os.path.exists(filepath):
            print(f"❌ File not found: {filepath}")
            continue
            
        try:
            parser = XMLParser(filepath)
            model = parser.parse()
            
            optimizer = DOFOptimizer(model)
            num_eq, semi_bw, full_bw = optimizer.optimize()
            
            print(f"{model.name:<35} | {num_eq:<12} | {semi_bw:<12} | {full_bw:<12}")
        except Exception as e:
            print(f"❌ Error processing {filepath}: {e}")

if __name__ == "__main__":
    run_tests()