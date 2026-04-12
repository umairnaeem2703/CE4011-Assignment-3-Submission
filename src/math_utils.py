# src/math_utils.py

def zeros(rows: int, cols: int) -> list:
    """Creates a 2D list of zeros."""
    return [[0.0 for _ in range(cols)] for _ in range(rows)]

def add(A: list, B: list) -> list:
    """Element-wise addition of two 2D matrices."""
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def subtract(A: list, B: list) -> list:
    """Element-wise subtraction of two 2D matrices."""
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def matmul(A: list, B: list) -> list:
    """Multiplies two 2D matrices (A * B)."""
    rows_A, cols_A = len(A), len(A[0])
    rows_B, cols_B = len(B), len(B[0])
    
    if cols_A != rows_B:
        raise ValueError("Incompatible matrix dimensions for multiplication.")
        
    result = zeros(rows_A, cols_B)
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                result[i][j] += A[i][k] * B[k][j]
    return result

def transpose(A: list) -> list:
    """Transposes a 2D matrix."""
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]

def invert_matrix(A: list) -> list:
    """Inverts a square 2D matrix using Gauss-Jordan elimination."""
    n = len(A)
    if n != len(A[0]):
        raise ValueError("Matrix must be square to invert.")
        
    # Create augmented matrix [A | I]
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(A)]
    
    for i in range(n):
        # Find pivot
        pivot_row = max(range(i, n), key=lambda r: abs(aug[r][i]))
        if abs(aug[pivot_row][i]) < 1e-12:
            raise ValueError("Matrix is singular or nearly singular.")
            
        # Swap rows
        aug[i], aug[pivot_row] = aug[pivot_row], aug[i]
        
        # Normalize pivot row
        pivot = aug[i][i]
        for j in range(i, 2 * n):
            aug[i][j] /= pivot
            
        # Eliminate other rows
        for k in range(n):
            if k != i:
                factor = aug[k][i]
                for j in range(i, 2 * n):
                    aug[k][j] -= factor * aug[i][j]
                    
    # Extract the inverted matrix (right half)
    return [row[n:] for row in aug]