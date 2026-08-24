import numpy as np
from functnTV import functnTV

def gradFunctnTV(P, betaTV, nprismsY):
    """
    Numerical Gradient of 2D TV.
    """
    tempP = np.zeros(len(P))
    grad = np.zeros(len(P))
    for k in range(len(P)):
        tempP[k] = P[k]
    
    for k in range(len(tempP)):
        aj = tempP[k]
        delta = 0.001 * tempP[k]
        if delta < 0.0001:
            delta = 0.001
        
        tempP[k] = aj + delta
        yf = functnTV(tempP, betaTV, nprismsY)
        
        tempP[k] = aj - delta
        td = functnTV(tempP, betaTV, nprismsY)
        
        grad[k] = (yf - td) / (2 * delta)
        tempP[k] = aj
        
    return grad
