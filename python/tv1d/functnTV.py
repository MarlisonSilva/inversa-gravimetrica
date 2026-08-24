import numpy as np


def functnTV(P, beta):
    """
    Computes the Total Variation functional of vector P.

    The functional is a smooth approximation of the L1-norm of the
    differences between adjacent elements, given by:

        TV(P) = sum_i sqrt((P[i+1] - P[i])^2 + beta^2)

    Parameters
    ----------
    P : np.ndarray
        Solution vector (prism thicknesses).
    beta : float
        Differentiability restoration parameter (small positive value).

    Returns
    -------
    TV : float
        Total variation value.
    """
    TV = 0.0
    for i in range(len(P) - 1):
        dp = P[i + 1] - P[i]
        TV += np.sqrt(dp ** 2 + beta ** 2)
    return TV
