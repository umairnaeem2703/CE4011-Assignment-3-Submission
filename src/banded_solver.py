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
        """
        Solves [K]{D} = {F} using in-place banded Gaussian elimination.
        Returns the global displacement vector {D}.
        """
        if self.num_eq == 0:
            return []

        # ==========================================
        # 1. FORWARD ELIMINATION
        # ==========================================
        for k in range(self.num_eq):
            pivot = self.K[k][0]
            
            # The Safety Net (Assignment Q3 requirement)
            # If the pivot is practically zero, the structure is a mechanism.
            if abs(pivot) < 1e-10:
                raise UnstableStructureError(
                    f"Mechanism detected! Instability found at Equation/DOF {k}. "
                    "Ensure your structure is fully restrained and adequately braced."
                )

            # Eliminate entries below the pivot (within the bandwidth)
            # Since the matrix is symmetric, K_ik = K_ki
            for i in range(k + 1, min(self.num_eq, k + self.semi_bw)):
                # The multiplier is K_ki / K_kk. 
                # In our upper-banded array, K_ki is stored at K[k][i - k]
                multiplier = self.K[k][i - k] / pivot
                
                # Update the load vector
                self.F[i][0] -= multiplier * self.F[k][0]
                
                # Update the remaining row entries in the stiffness matrix
                for j in range(i, min(self.num_eq, k + self.semi_bw)):
                    # K_ij = K_ij - multiplier * K_kj
                    self.K[i][j - i] -= multiplier * self.K[k][j - k]

        # ==========================================
        # 2. BACK SUBSTITUTION
        # ==========================================
        # Initialize displacement vector {D}
        D = [[0.0] for _ in range(self.num_eq)]
        
        for i in range(self.num_eq - 1, -1, -1):
            sum_val = self.F[i][0]
            
            # Substitute known displacements from the upper band
            for j in range(i + 1, min(self.num_eq, i + self.semi_bw)):
                sum_val -= self.K[i][j - i] * D[j][0]
                
            D[i][0] = sum_val / self.K[i][0]

        return D