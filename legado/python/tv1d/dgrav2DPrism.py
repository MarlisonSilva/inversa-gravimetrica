import numpy as np

def dgrav2DPrism(den, x0, z0, w, h, xobs, unit):
    """
    Derivada Analítica da Anomalia em relação à espessura h.
    
    Parâmetros: den, x0, z0, w, h, xobs, unit (0=km, 1=m)
    """
    CG = 6.672e-08
    cd = den * 1000.0
    x = xobs - x0
    
    p1 = 2 * (z0 + h) * (1.0 / ((z0 + h)**2 + x**2) - 1.0 / ((z0 + h)**2 + (x - w)**2))
    p2 = 2 * (z0 + h) / ((z0 + h)**2 + (x - w)**2)
    
    t1 = np.arctan2(x - w, z0 + h)
    t2 = np.arctan2(x, z0 + h)
    
    p3 = t1 - t2 + x / ((z0 + h) * (1.0 + x**2 / (z0 + h)**2)) - \
         (x - w) / ((z0 + h) * (1.0 + (x - w)**2 / (z0 + h)**2))
         
    dg = 2 * 10**5 * CG * cd * ((x / 2.0) * p1 + (w / 2.0) * p2 - p3)
    
    if unit == 1:
        dg = dg / 1000.0
        
    return dg
