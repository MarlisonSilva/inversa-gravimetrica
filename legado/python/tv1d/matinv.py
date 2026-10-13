import numpy as np

def matinv(ARRAY, NORDER):
    """
    Inversão de Matriz por eliminação de Gauss-Jordan.
    Convertido fielmente do MATLAB. 
    NOTA: Em Python moderno, seria preferível usar `np.linalg.inv(ARRAY)`.
    """
    A = np.copy(ARRAY)
    SFLAG = 0
    DET = 1.0
    IK = np.zeros(NORDER, dtype=int)
    JK = np.zeros(NORDER, dtype=int)
    
    for K in range(NORDER):
        AMAX = 0.0
        for I in range(K, NORDER):
            for J in range(K, NORDER):
                if (abs(AMAX) - abs(A[I, J])) <= 0:
                    AMAX = A[I, J]
                    IK[K] = I
                    JK[K] = J
                    
        if AMAX == 0.0:
            DET = 0.0
            SFLAG = 1
            return A, SFLAG
            
        I = IK[K]
        if (I - K) > 0:
            for J in range(NORDER):
                SAVE1 = A[K, J]
                A[K, J] = A[I, J]
                A[I, J] = -SAVE1
                
        J = JK[K]
        if (J - K) > 0:
            for I_idx in range(NORDER):
                SAVE1 = A[I_idx, K]
                A[I_idx, K] = A[I_idx, J]
                A[I_idx, J] = -SAVE1
                
        for I_idx in range(NORDER):
            if (I_idx - K) != 0:
                A[I_idx, K] = -A[I_idx, K] / AMAX
                
        for I_idx in range(NORDER):
            for J_idx in range(NORDER):
                if (I_idx - K) != 0 and (J_idx - K) != 0:
                    A[I_idx, J_idx] = A[I_idx, J_idx] + A[I_idx, K] * A[K, J_idx]
                    
        for J_idx in range(NORDER):
            if (J_idx - K) != 0:
                A[K, J_idx] = A[K, J_idx] / AMAX
                
        A[K, K] = 1.0 / AMAX
        
    # Reordering loop
    for L in range(NORDER):
        K = NORDER - 1 - L
        J = IK[K]
        if (J - K) > 0:
            for I_idx in range(NORDER):
                SAVE1 = A[I_idx, K]
                A[I_idx, K] = -A[I_idx, J]
                A[I_idx, J] = SAVE1
                
        I = JK[K]
        if (I - K) > 0:
            for J_idx in range(NORDER):
                SAVE1 = A[K, J_idx]
                A[K, J_idx] = -A[I, J_idx]
                A[I, J_idx] = SAVE1
                
    return A, SFLAG
