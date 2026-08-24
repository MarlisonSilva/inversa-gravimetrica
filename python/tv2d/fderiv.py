import numpy as np
from gpDeriv import gpDeriv

def fderiv(i, j, P, anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprismsY, prismDensity):
    """
    Partial Derivative 2D (Jacobian row (i,j)).
    i, j must be 0-based in Python.
    """
    deriv = np.zeros(len(P))
    w = j + i * nobsY
    
    for k in range(len(P)):
        row = int(k // nprismsY)
        col = int(k % nprismsY)
        xcr = (row * prismWidthX + (row + 1) * prismWidthX) / 2
        ycr = (col * prismWidthY + (col + 1) * prismWidthY) / 2
        
        deriv[k] = gpDeriv(prismDensity, anomObs[w, 0], anomObs[w, 1], 0.0001 + P[k], 0.0001, xcr, ycr, prismWidthX/2, prismWidthY/2)
        
    return deriv
