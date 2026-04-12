from parser import XMLParser
from dof_optimizer import DOFOptimizer
from matrix_assembly import MatrixAssembler
from banded_solver import BandedSolver, UnstableStructureError

# Load a structurally sound model
parser = XMLParser("./data/example1_case1_truss.xml")
model = parser.parse()

opt = DOFOptimizer(model)
num_eq, semi_bw, _ = opt.optimize()

assembler = MatrixAssembler(model, num_eq, semi_bw)
K_banded, F_global = assembler.assemble("LC1")

solver = BandedSolver(K_banded, F_global, semi_bw)

try:
    displacements = solver.solve()
    print("✅ System solved successfully!")
    print("Displacements:", [d[0] for d in displacements[:]])
except UnstableStructureError as e:
    print(f"❌ {e}")