import numpy as np

def createInitialPrisms(nobs, P, udist, anomObs, prismWidth, prismDensity):
    """
    Modelo Inicial por Aproximação de Bouguer.
    """
    ufact = 1e3
    if udist == 0:
        ufact = 1.0
        
    dx = prismWidth
    rzero = dx / 4.0
    
    for k in range(len(P)):
        xcp = k * dx + dx / 2.0
        P[k] = 0.0
        r = rzero
        
        while True:
            m = 0
            anom = 0.0
            weight = 0.0
            
            for i in range(nobs):
                d = abs(xcp - anomObs[i, 0])
                if (d < r) and (anomObs[i, 1] != -9999.00):
                    anom = anom + anomObs[i, 1] / (d**2 + 0.0001)
                    m = m + 1
                    weight = weight + 1.0 / (d**2 + 0.0001)
                    
            if m != 0:
                anom = anom / weight
                P[k] = ufact * anom / (41.9 * prismDensity + 0.0001)
                break
            else:
                r = r + rzero
                
    return P
