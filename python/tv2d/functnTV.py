import numpy as np

def suavy(n, nx):
    """
    Creates an index array for reordering prisms by column (Y).
    MATLAB uses 1-based, Python uses 0-based indexing.
    """
    m = np.zeros(n, dtype=int)
    k = 0
    for j in range(nx):
        for i in range(j, n, nx):
            m[k] = i
            k += 1
    return m

def functnTV(P, betaTV, nprismsY):
    """
    2D Total Variation Functional.
    """
    m = suavy(len(P), nprismsY)
    fTV = 0.0
    for i in range(len(P) - 1):
        fTV += np.sqrt(np.abs(P[i+1] - P[i])**2 + betaTV**2)
        fTV += np.sqrt(np.abs(P[m[i+1]] - P[m[i]])**2 + betaTV**2)
    return fTV
