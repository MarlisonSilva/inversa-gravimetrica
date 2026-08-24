import numpy as np
try:
    from .dgrav2DPrism import dgrav2DPrism
except ImportError:
    from dgrav2DPrism import dgrav2DPrism

def fderiv(i, P, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, unit):
    """
    Derivada Parcial do Funcional (Jacobiana linha i).
    i é 0-based index no Python para referenciar XOBS[i].
    """
    DERIV = np.zeros(len(P))
    
    if dcType == 0:
        for k in range(len(P)):
            xc = xminPrisms + k * prismWidth + prismWidth / 2.0
            DERIV[k] = dgrav2DPrism(prismDensity, xc, 0.0005, prismWidth, P[k], XOBS[i], unit)
    elif dcType == 1:
        raise NotImplementedError("dcType=1 (hipérbole) não está implementado.")
        
    return DERIV
