import numpy as np
import scipy.io
import os

try:
    from .createInitialPrisms import createInitialPrisms
    from .functn import functn
    from .QNiterate import QNiterate
except ImportError:
    from createInitialPrisms import createInitialPrisms
    from functn import functn
    from QNiterate import QNiterate

def inv_2d(anomObs, maxIter=100, flambda=0.001, usrTOL=0.0001,
           alpha=1.25, betaTV=0.001,
           nprisms_x=42, nprisms_y=26,
           xmin=0, xmax=52.5, ymin=0, ymax=32.5,
           prismDensity=-0.30, dcType=0, dcBeta=0.0,
           nobsX=21, nobsY=13):
    """
    Complete 2D Inversion.
    """
    n_params = nprisms_x * nprisms_y
    P = np.zeros(n_params)
    
    print("Creating initial model...")
    P = createInitialPrisms(P, prismDensity, anomObs, nobsX, nobsY, nprisms_x, nprisms_y, xmax, xmin, ymax, ymin)
    
    prismWidthX = (xmax - xmin) / nprisms_x
    prismWidthY = (ymax - ymin) / nprisms_y
    
    YFIT = functn(P, anomObs, nobsX, nobsY, xmin, xmax, ymin, ymax, nprisms_x, nprisms_y, prismDensity)
    
    state = {'Hk': None, 'bbeta': None}
    
    print("Starting inversion...")
    for it in range(maxIter):
        chisqr, chisqA, YFIT, P, flambda, state = QNiterate(
            P, YFIT, anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprisms_y,
            prismDensity, dcType, dcBeta, alpha, betaTV, flambda, it,
            xmin, xmax, ymin, ymax, nprisms_x, nprisms_y, state
        )
        
        print(f"Iter: {it+1:3d} | Obj: {chisqr:.6f} | RMSE: {chisqA:.6f} | lambda: {flambda:.4f}")
        
        if chisqA < usrTOL:
            print("Converged by usrTOL criterion.")
            break
            
    return P

if __name__ == '__main__':
    mat_file = 'anomObs.mat'
    if os.path.exists(mat_file):
        mat = scipy.io.loadmat(mat_file)
        anomObs = mat['anomObs']
        inv_2d(anomObs)
    else:
        print(f"{mat_file} not found. Please provide the data to run the script.")
