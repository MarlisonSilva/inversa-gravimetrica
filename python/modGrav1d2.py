from grav2DPrism2 import grav2DPrism2

def modGrav1d2(P, XOBS, iobs, xminPrisms, prismWidth, prismDensity):
    """
    Computes the 2D gravity anomaly for an array of prisms.
    """
    functn = 0.0
    ztop = 0.0005
    
    for i in range(len(P)):
        # xc is the center of the i-th prism
        xc = xminPrisms + i * prismWidth + prismWidth / 2.0
        # P[i] is the height/thickness of the i-th prism
        functn += grav2DPrism2(prismDensity, xc, ztop, prismWidth, P[i], XOBS[iobs])
        
    functn = functn / 1000.0
    return functn
