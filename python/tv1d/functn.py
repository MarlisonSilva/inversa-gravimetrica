import numpy as np
try:
    from .grav2DPrism import grav2DPrism
except ImportError:
    from grav2DPrism import grav2DPrism

def functn(P, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, unit):
    """
    Forward Gravity (modelo → anomalia) 1D
    """
    ztop = 0.0005
    misfit = np.zeros(len(XOBS))
    
    for i in range(len(P)):
        xc = xminPrisms + i * prismWidth + prismWidth / 2.0
        for j in range(len(XOBS)):
            if dcType == 0:
                misfit[j] += grav2DPrism(prismDensity, xc, ztop, prismWidth, P[i], XOBS[j])
            elif dcType == 1:
                raise NotImplementedError("dcType=1 (hipérbole) não está implementado.")
                
    if unit == 1:
        misfit = misfit / 1000.0
        
    return misfit
