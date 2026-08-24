import numpy as np
try:
    from .functnTV import functnTV
except ImportError:
    from functnTV import functnTV


def fchisq(GOBS, YFIT, P, alpha, betaTV):
    """
    Computes the objective function for the TV-regularized inversion:

        T(P) = ||Y_obs - F(P)||^2 + alpha * TV(P)

    Parameters
    ----------
    GOBS : np.ndarray
        Observed gravity data (shape: nobs).
    YFIT : np.ndarray
        Gravity data computed from current model P (shape: nobs).
    P : np.ndarray
        Current solution vector (prism thicknesses).
    alpha : float
        Regularization parameter for the Total Variation functional.
    betaTV : float
        Differentiability restoration parameter for TV.

    Returns
    -------
    chisqR : float
        Total objective function value (misfit + alpha * TV).
    chisqA : float
        Normalized RMS misfit (sqrt of mean squared residual).
    tvValue : float
        Total Variation value of current model.
    """
    npts = len(GOBS)
    chisq = np.sum((GOBS - YFIT) ** 2)

    chisqA = np.sqrt(chisq / npts)
    tvValue = functnTV(P, betaTV)
    chisqR = chisq + alpha * tvValue

    return chisqR, chisqA, tvValue
