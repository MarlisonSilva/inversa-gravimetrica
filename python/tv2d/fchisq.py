import numpy as np
from functnTV import functnTV

def fchisq(anomObs, YFIT, P, alpha, betaTV, nprismsY, nobsX, nobsY):
    """
    2D Objective Function.
    """
    chisq = 0.0
    k = 0
    
    for i in range(nobsX):
        for j in range(nobsY):
            w = j + i * nobsY
            if anomObs[w, 2] != -9999.00:
                chisq += (anomObs[w, 2] - YFIT[i, j])**2
                k += 1
                
    if k > 0:
        chisqA = np.sqrt(chisq / k)
    else:
        chisqA = 0.0
        
    Obj = chisq + alpha * functnTV(P, betaTV, nprismsY)
    
    return Obj, chisqA
