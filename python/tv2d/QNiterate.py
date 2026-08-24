import numpy as np
from fderiv import fderiv
from gradFunctnTV import gradFunctnTV
from fchisq import fchisq
from functn import functn
from chodec import chodec

def QNiterate(P, YFIT, anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprismsY,
              prismDensity, dcType, dcBeta, alpha, betaTV, flambda, it,
              xmin, xmax, ymin, ymax, nprisms_x, nprisms_y, state):
    """
    Quasi-Newton Iteration 2D with TV.
    state: dict with 'Hk' and 'bbeta', initialized in inv_2d
    """
    n_params = len(P)
    I = np.eye(n_params)
    
    # Init state if first iteration
    if state.get('Hk') is None:
        state['Hk'] = I.copy()
    if state.get('bbeta') is None:
        state['bbeta'] = I.copy()
        
    Hk = state['Hk']
    bbeta = state['bbeta']
    
    # Calculate gradient
    grad = gradFunctnTV(P, betaTV, nprismsY)
    
    # Construct Jacobian A
    A = np.zeros((nobsX * nobsY, n_params))
    for i in range(nobsX):
        for j in range(nobsY):
            w = j + i * nobsY
            if anomObs[w, 2] != -9999.00:
                deriv = fderiv(i, j, P, anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprismsY, prismDensity)
                A[w, :] = deriv
                
    dy = np.zeros(nobsX * nobsY)
    for i in range(nobsX):
        for j in range(nobsY):
            w = j + i * nobsY
            if anomObs[w, 2] != -9999.00:
                dy[w] = anomObs[w, 2] - YFIT[i, j]
                
    gradObj = -A.T @ dy + alpha * grad
    
    # Marquardt term
    H = A.T @ A + alpha * bbeta
    diagH = np.diag(np.diag(H))
    
    H_lambda = H + flambda * diagH
    
    invH, sflag = chodec(H_lambda)
    if sflag != 0:
        # fallback to standard inversion if cholesky fails
        invH = np.linalg.inv(H_lambda)
        
    deltaP = -invH @ gradObj
    
    # Try step
    newP = P + deltaP
    newYFIT = functn(newP, anomObs, nobsX, nobsY, xmin, xmax, ymin, ymax, nprisms_x, nprisms_y, prismDensity)
    newObj, newChisqA = fchisq(anomObs, newYFIT, newP, alpha, betaTV, nprismsY, nobsX, nobsY)
    
    Obj, chisqA = fchisq(anomObs, YFIT, P, alpha, betaTV, nprismsY, nobsX, nobsY)
    
    if newObj < Obj:
        P = newP
        YFIT = newYFIT
        flambda = flambda / 10.0
        chisqr = newObj
        chisqA = newChisqA
        
        # update bbeta (Hessian approximation BFGS logic)
        sk = deltaP
        yk = (A.T @ A + alpha * bbeta) @ sk
        if np.dot(yk, sk) > 0:
            rho = 1.0 / np.dot(yk, sk)
            V = I - rho * np.outer(sk, yk)
            Hk = V.T @ Hk @ V + rho * np.outer(sk, sk)
            bbeta = np.linalg.inv(Hk) # keeping approximation
            
    else:
        flambda = flambda * 10.0
        chisqr = Obj
        
    state['Hk'] = Hk
    state['bbeta'] = bbeta
    
    return chisqr, chisqA, YFIT, P, flambda, state
