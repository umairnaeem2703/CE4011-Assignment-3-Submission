# src/banded_solver.py

class UnstableStructureError(Exception):
    """Custom exception raised when a structure is a mechanism."""
    pass

class BandedSolver:
    def __init__(self, K_banded: list, F_global: list, semi_bandwidth: int):
        self.K = K_banded
        self.F = F_global
        self.semi_bw = semi_bandwidth
        self.num_eq = len(K_banded)

    def solve(self) -> list:
        """Solves [K]{D} = {F} using in-place banded Gaussian elimination."""
        if self.num_eq == 0:
            return []

        F = [[self.F[i][0]] for i in range(self.num_eq)]
        K = [row[:] for row in self.K]

        # Forward elimination
        for k in range(self.num_eq):
            pivot = K[k][0]
            
            if abs(pivot) < 1e-10:
                raise UnstableStructureError(
                    f"Mechanism detected at Equation/DOF {k}. "
                    "Ensure structure is fully restrained and adequately braced."
                )

            for i in range(k + 1, min(self.num_eq, k + self.semi_bw)):
                multiplier = K[k][i - k] / pivot
                F[i][0] -= multiplier * F[k][0]
                
                for j in range(i, min(self.num_eq, k + self.semi_bw)):
                    K[i][j - i] -= multiplier * K[k][j - k]

        # Back substitution
        D = [[0.0] for _ in range(self.num_eq)]
        
        for i in range(self.num_eq - 1, -1, -1):
            sum_val = F[i][0]
            for j in range(i + 1, min(self.num_eq, i + self.semi_bw)):
                sum_val -= K[i][j - i] * D[j][0]
            D[i][0] = sum_val / K[i][0]

        return D