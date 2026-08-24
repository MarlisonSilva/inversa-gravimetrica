import numpy as np

try:
    from .functnTV import functnTV
    from .gradFunctnTV import gradFunctnTV
    from .fchisq import fchisq
    from .functn import functn
    from .fderiv import fderiv
    from .matinv import matinv
except ImportError:
    from functnTV import functnTV
    from gradFunctnTV import gradFunctnTV
    from fchisq import fchisq
    from functn import functn
    from fderiv import fderiv
    from matinv import matinv


def QNiterate(nprisms, P, anomObs, YFIT, xminPrisms, prismWidth, prismDensity,
              dcType, dcBeta, alpha, betaTV, flambda, it, unit, state):
    """
    Quasi-Newton iteration (BFGS) with Total Variation regularization for 1D gravity inversion.

    The MATLAB 'persistent' variables Hk and bbeta are replaced by the `state` dict,
    which must be initialized externally before the first call:
        state = {'Hk': np.eye(nprisms), 'bbeta': np.zeros(nprisms)}

    Parameters
    ----------
    nprisms : int
        Number of prisms.
    P : np.ndarray
        Current model (prism thicknesses).
    anomObs : np.ndarray
        Observed data, shape (nobs, 2) — columns: [x, gz].
    YFIT : np.ndarray
        Current forward response.
    xminPrisms : float
        Start of the discretization along the profile.
    prismWidth : float
        Horizontal width of each prism.
    prismDensity : float
        Density contrast of prisms (g/cc).
    dcType : int
        Density contrast type (0=constant, 1=hyperbolic).
    dcBeta : float
        Decaying factor (used if dcType=1).
    alpha : float
        TV regularization parameter.
    betaTV : float
        TV differentiability restoration parameter.
    flambda : float
        Marquardt (damping) parameter.
    it : int
        Current iteration number (1-based, matching MATLAB).
    unit : int
        Distance unit (0=km, 1=m).
    state : dict
        Mutable state dict with keys 'Hk' (np.ndarray) and 'bbeta' (np.ndarray).

    Returns
    -------
    chisqr : float
    chisqA : float
    tvValue : float
    YFIT : np.ndarray
    P : np.ndarray
    flambda : float
    state : dict
    """
    nterms = len(P)
    Hk = state['Hk']
    bbeta = state['bbeta']

    XOBS = anomObs[:, 0]
    GOBS = anomObs[:, 1]
    nobs = len(XOBS)

    B = np.zeros(nterms)
    idfp = 0  # 0 = use BFGS, 1 = use DFP

    # --- Compute gradient ---
    if it == 1:
        beta1 = np.zeros(nterms)
        for i in range(nobs):
            DERIV = fderiv(i, P, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, unit)
            beta1 += (GOBS[i] - YFIT[i]) * DERIV
        beta3 = gradFunctnTV(P, betaTV)
        beta = -2.0 * beta1 + alpha * beta3
    else:
        beta = bbeta.copy()

    chisq1, chisqA, tvValue = fchisq(GOBS, YFIT, P, alpha, betaTV)

    # --- Marquardt loop ---
    lambda_ = flambda
    while True:
        ARRAY = np.zeros((nterms, nterms))
        for j in range(nterms):
            for k in range(nterms):
                denom = np.sqrt(abs(Hk[j, j] * Hk[k, k]))
                ARRAY[j, k] = Hk[j, k] / denom if denom != 0 else 0.0
            ARRAY[j, j] += lambda_

        # NOTE: np.linalg.inv(ARRAY) is the modern equivalent.
        ARRAY, SFlag = matinv(ARRAY, nterms)

        for j in range(nterms):
            B[j] = P[j]
            for k in range(nterms):
                denom = np.sqrt(abs(Hk[j, j] * Hk[k, k]))
                B[j] -= beta[k] * ARRAY[j, k] / denom if denom != 0 else 0.0
            if B[j] < 0:
                B[j] = 0.01

        YFIT_new = functn(B, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, unit)
        chisqr, chisqA, tvValue = fchisq(GOBS, YFIT_new, B, alpha, betaTV)

        if (chisq1 - chisqr) < 0:
            if lambda_ > 1e8:
                flambda = -1.0
                state['Hk'] = Hk
                state['bbeta'] = bbeta
                return chisqr, chisqA, tvValue, YFIT_new, P, flambda, state
            lambda_ *= 10.0
            continue
        else:
            # Accept solution — compute new gradient
            beta1_new = np.zeros(nterms)
            for i in range(nobs):
                DERIV = fderiv(i, B, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, unit)
                beta1_new += (GOBS[i] - YFIT_new[i]) * DERIV
            bbeta3 = gradFunctnTV(B, betaTV)
            bbeta = -2.0 * beta1_new + alpha * bbeta3

            SK = B - P
            YK = bbeta - beta

            # BFGS Hessian update
            XM2 = np.outer(YK, YK)
            E1 = np.dot(YK, SK)
            V1 = Hk @ SK
            XM1 = np.outer(V1, V1)
            E2 = np.dot(SK, V1)

            eps = 1e-6
            Hk = Hk + XM2 / (E1 + eps) - XM1 / (E2 + eps)

            P = B.copy()
            YFIT = YFIT_new
            state['Hk'] = Hk
            state['bbeta'] = bbeta
            return chisqr, chisqA, tvValue, YFIT, P, flambda, state
