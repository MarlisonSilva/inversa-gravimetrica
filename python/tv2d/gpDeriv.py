import numpy as np
from f_dkern import f_dkern

def gpDeriv(rho, xr, yr, hb, ht, xc, yc, xw2, yw2):
    """
    Analytical Derivative of the 3D Prism relative to thickness.
    """
    GAMMA = 6.67
    x1 = xc - xr - xw2
    x2 = xc - xr + xw2
    y1 = yc - yr - yw2
    y2 = yc - yr + yw2
    z1 = ht
    z2 = hb
    if z2 == 0:
        z2 = 1e3
        
    pDeriv = GAMMA * rho * (f_dkern(x1, y2, z2) + f_dkern(x2, y1, z2) -
                            f_dkern(x2, y2, z2) - f_dkern(x1, y1, z2))
    return pDeriv
