import numpy as np

def createInitialPrisms(P, prismDensity, anomObs, nobsX, nobsY, nprisms_x, nprisms_y, xmax, xmin, ymax, ymin):
    """
    2D Initial Model.
    """
    EPS = 1e-4
    dx = (xmax - xmin) / nprisms_x
    dy = (ymax - ymin) / nprisms_y
    rzero = dx / 4
    
    for k in range(len(P)):
        row = int(k // nprisms_y)
        col = int(k % nprisms_y)
        xcp = (row * dx + (row + 1) * dx) / 2
        ycp = (col * dy + (col + 1) * dy) / 2
        P[k] = 0.0
        
        r = rzero
        while True:
            m = 0
            campo = 0.0
            peso = 0.0
            for i in range(nobsX):
                for j in range(nobsY):
                    w = j + i * nobsY
                    d = np.sqrt((xcp - anomObs[w, 0])**2 + (ycp - anomObs[w, 1])**2)
                    if d < r and anomObs[w, 2] != -9999.00:
                        campo += anomObs[w, 2] / (d**2 + EPS)
                        m += 1
                        peso += 1.0 / (d**2 + EPS)
                        
            if m == 0:
                r += rzero
            else:
                campo = campo / peso
                P[k] = campo / (41.9 * prismDensity + EPS)
                break
                
    return P
