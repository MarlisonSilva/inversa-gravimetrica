import numpy as np

def grav2DPrism(den, x0, z0, w, h, xobs):
    """
    Anomalia Gravimétrica Analítica 2D (caso 1D, prisma isolado).
    """
    CG = 6.672e-08
    cd = den * 1000.0
    x = xobs - x0 + w / 2.0
    
    r1 = np.sqrt(x**2 + z0**2)
    r2 = np.sqrt(x**2 + (h + z0)**2)
    r3 = np.sqrt((x - w)**2 + z0**2)
    r4 = np.sqrt((x - w)**2 + (h + z0)**2)
    
    # Prevenindo divisão por zero usando arctan2 se z0 == 0
    theta1 = np.arctan2(x, z0) if z0 != 0 else (np.pi/2 * np.sign(x) if x != 0 else 0)
    theta2 = np.arctan2(x, z0 + h)
    theta3 = np.arctan2(x - w, z0) if z0 != 0 else (np.pi/2 * np.sign(x - w) if (x-w) != 0 else 0)
    theta4 = np.arctan2(x - w, z0 + h)
    
    gz = 2 * 10**5 * CG * cd * (
        x * np.log((r2 * r3) / (r1 * r4)) + 
        w * np.log(r4 / r3) - 
        h * (theta4 - theta2) + 
        z0 * (theta3 - theta1)
    )
    return gz
