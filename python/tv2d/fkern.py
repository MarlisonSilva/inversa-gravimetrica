import numpy as np

def FKERN(X, Y, Z):
    """
    Auxiliary Function FKERN (potential kernel).
    """
    R = np.sqrt(X*X + Y*Y + Z*Z)
    CREAL = Z
    CIMAG = R + X + Y
    PHI = np.arctan2(CIMAG, CREAL)
    
    if X == 0:
        Kern = Y * np.log(R + X) + 2 * Z * PHI
    elif Y == 0:
        Kern = X * np.log(R + Y) + 2 * Z * PHI
    elif X == 0 and Y == 0:
        Kern = 2 * Z * PHI
    else:
        Kern = X * np.log(R + Y) + Y * np.log(R + X) + 2 * Z * PHI
        
    return Kern
