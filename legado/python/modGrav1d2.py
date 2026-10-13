from grav2DPrism2 import grav2DPrism2

def modGrav1d2(P, XOBS, iobs, xminPrisms, prismWidth, prismDensity, unit=0):
    """
    Computes the 2D gravity anomaly for an array of prisms.

    unit: 0 - distances in kilometers, result in mGal.
          1 - distances in meters, result converted from microGal to mGal.
    """
    functn = 0.0
    ztop = 0.0005

    for i in range(len(P)):
        # xc is the center of the i-th prism
        xc = xminPrisms + i * prismWidth + prismWidth / 2.0
        # P[i] is the height/thickness of the i-th prism
        functn += grav2DPrism2(prismDensity, xc, ztop, prismWidth, P[i], XOBS[iobs])

    # Conversion from microGal to mGal if distances are in meters.
    if unit == 1:
        functn = functn / 1000.0

    return functn
