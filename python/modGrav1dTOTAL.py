from modGrav1d2 import modGrav1d2

def modGrav1dTOTAL(P, XOBS, xminPrisms, prismWidth, prismDensity):
    """
    Computes the total 1D gravity anomaly for a series of observation points.
    """
    Fmod1 = [0.0] * len(XOBS)
    
    for iobs in range(len(XOBS)):
        Fmod1[iobs] = modGrav1d2(P, XOBS, iobs, xminPrisms, prismWidth, prismDensity)
        
    return Fmod1
