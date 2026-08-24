import numpy as np

def f_dkern(x, y, z):
    """
    Derivative of the Potential Kernel.
    """
    R = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(R + x + y, z)
    
    # 2*z*(...) term expression converted faithfully
    term = 2 * z * ((1/R - (1/z**2)*(R+x+y)) / (1 + (1/z**2)*(R+x+y)**2))
    
    if x == 0:
        dkern = y * z / (R * (R + x)) + 2 * phi + term
    elif y == 0:
        dkern = x * z / (R * (R + y)) + 2 * phi + term
    elif x == 0 and y == 0:
        dkern = 2 * phi + term
    else:
        dkern = x * z / (R * (R + y)) + y * z / (R * (R + x)) + 2 * phi + term
        
    return dkern
