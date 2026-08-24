import numpy as np

def chodec(A):
    """
    Inversion via Cholesky Decomposition.
    This is a faithful manual translation of the original CHODEC.m.
    Note: np.linalg.inv(A) or np.linalg.cholesky(A) are the modern Python alternatives.
    Returns (A_inv, SFLAG) where SFLAG=0 if OK.
    """
    n = A.shape[0]
    SFLAG = 0
    IFLAG = 0
    
    # 1. if A[0,0] < 0, multiply by -1
    if A[0, 0] < 0:
        IFLAG = 1
        A = -A
        
    # 2. Lower triangular decomposition
    L = np.zeros_like(A, dtype=float)
    for i in range(n):
        for j in range(i + 1):
            s = sum(L[i, k] * L[j, k] for k in range(j))
            if i == j:
                val = A[i, i] - s
                if val <= 0:
                    SFLAG = 1
                    return A, SFLAG
                L[i, j] = np.sqrt(val)
            else:
                L[i, j] = (1.0 / L[j, j]) * (A[i, j] - s)
                
    # 3. Compute inverse by substitution
    invL = np.zeros_like(L, dtype=float)
    for i in range(n):
        invL[i, i] = 1.0 / L[i, i]
        for j in range(i):
            s = sum(L[i, k] * invL[k, j] for k in range(j, i))
            invL[i, j] = -s / L[i, i]
            
    # 4. Symmetrize
    A_inv = invL.T @ invL
    
    # 5. Undo sign if IFLAG
    if IFLAG == 1:
        A_inv = -A_inv
        
    return A_inv, SFLAG
