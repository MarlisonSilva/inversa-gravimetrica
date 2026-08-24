import numpy as np
from gprism import GPRISM

def functn(P, anomObs, nobsX, nobsY, xmin, xmax, ymin, ymax, nprisms_x, nprisms_y, prismDensity):
    """
    Forward Gravity 2D (model -> anomaly).
    """
    htop = 0.0001
    dx = (xmax - xmin) / nprisms_x
    dy = (ymax - ymin) / nprisms_y
    Grav = np.zeros((nobsX, nobsY))
    
    for i in range(nobsX):
        for j in range(nobsY):
            # w in Python logic 0-based
            w = j + i * nobsY
            xob = anomObs[w, 0]
            yob = anomObs[w, 1]
            if anomObs[w, 2] != -9999.00:
                for k in range(len(P)):
                    row = int(k // nprisms_y)
                    col = int(k % nprisms_y)
                    hbase = P[k] + htop
                    xcr = (row * dx + (row + 1) * dx) / 2
                    ycr = (col * dy + (col + 1) * dy) / 2
                    Grav[i, j] += GPRISM(xob, yob, dx/2, dy/2, htop, hbase, xcr, ycr, prismDensity)
            else:
                Grav[i, j] = -9999.00
                
    return Grav
