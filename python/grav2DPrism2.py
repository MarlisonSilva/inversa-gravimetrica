import math

def grav2DPrism2(den, x0, z0, w, h, xobs):
    """
    Computes the gravity field for a 2D prism with constant density contrast.
    Ref.: Telford, Geldart, Sheriff. Applied Geophysics p. 73, eq. 2.49c.
    
    Parameters:
      den  Input - Density contrast (g/cc)
      x0   Input - Horizontal centre of the prism
      z0   Input - Top of the prism
      w    Input - Prism width
      h    Input - Prism thickness
      xobs Input - Observational point where to compute the field
      
    Return:
      Value of the field (mGal or microGal depending on inputs)
    """
    CG = 6.672e-08
    cd = den * 1000.0

    x = xobs - x0 + w / 2.0

    r1 = math.sqrt(x**2 + z0**2)
    r2 = math.sqrt(x**2 + (h + z0)**2)
    r3 = math.sqrt((x - w)**2 + z0**2)
    r4 = math.sqrt((x - w)**2 + (h + z0)**2)
    
    theta1 = math.atan(x / z0) if z0 != 0 else math.pi / 2 * (1 if x > 0 else -1 if x < 0 else 0)
    theta2 = math.atan(x / (z0 + h))
    theta3 = math.atan((x - w) / z0) if z0 != 0 else math.pi / 2 * (1 if (x - w) > 0 else -1 if (x - w) < 0 else 0)
    theta4 = math.atan((x - w) / (z0 + h))

    gz = 2 * (10**5) * CG * cd * (x * math.log((r2 * r3) / (r1 * r4)) + w * math.log(r4 / r3) - h * (theta4 - theta2) + z0 * (theta3 - theta1))

    return gz
