import numpy as np
try:
    from .functnTV import functnTV
except ImportError:
    from functnTV import functnTV

def gradFunctnTV(P, beta):
    """
    Gradiente Numérico da Variação Total 1D.
    """
    grad = np.zeros(len(P))
    tempP = P.copy()
    for i in range(len(P)):
        aj = tempP[i]
        dp = 0.01 * tempP[i]
        if dp < 0.001:
            dp = 0.001
        tempP[i] = aj + dp
        a = functnTV(tempP, beta)
        tempP[i] = aj - dp
        b = functnTV(tempP, beta)
        grad[i] = (a - b) / (2 * dp)
        tempP[i] = aj
    return grad
