import numpy as np

try:
    from .createInitialPrisms import createInitialPrisms
    from .functn import functn
    from .fchisq import fchisq
    from .QNiterate import QNiterate
except ImportError:
    from createInitialPrisms import createInitialPrisms
    from functn import functn
    from fchisq import fchisq
    from QNiterate import QNiterate


def inv_1d(anomObs, maxIter=400, flambda=0.001, usrTOL=0.0001,
           alpha=1.0, betaTV=0.001,
           nprisms=100, xminPrisms=0.0, xmaxPrisms=200.0,
           prismDensity=-0.3, dcType=0, dcBeta=0.0, udist=0,
           noise_sigma=0.0):
    """
    1D Gravity inversion with Total Variation regularization.

    Uses a Quasi-Newton (BFGS) optimizer with Marquardt damping.

    Parameters
    ----------
    anomObs : np.ndarray, shape (N, 2)
        Observed data — columns: [x_position, gz_value].
        Use gz = -9999.0 to flag missing observations.
    maxIter : int
        Maximum number of iterations.
    flambda : float
        Initial Marquardt damping parameter.
    usrTOL : float
        Convergence tolerance on relative change of objective function.
    alpha : float
        Total Variation regularization parameter.
    betaTV : float
        TV differentiability restoration parameter.
    nprisms : int
        Number of prisms in the interpretation model.
    xminPrisms : float
        Start of discretization along profile.
    xmaxPrisms : float
        End of discretization along profile.
    prismDensity : float
        Density contrast of all prisms (g/cc).
    dcType : int
        Density contrast type (0=constant, 1=hyperbolic).
    dcBeta : float
        Decaying factor (used if dcType=1).
    udist : int
        Distance unit (0=km, 1=m).
    noise_sigma : float
        Standard deviation of Gaussian noise to add to observed data (0 = no noise).

    Returns
    -------
    P : np.ndarray
        Final model — prism thicknesses.
    """
    prismWidth = (xmaxPrisms - xminPrisms) / nprisms
    nobs = len(anomObs)

    # Add optional Gaussian noise
    if noise_sigma > 0.0:
        noise = noise_sigma * np.random.randn(nobs)
        anomObs = anomObs.copy()
        anomObs[:, 1] = anomObs[:, 1] + anomObs[:, 1] * noise

    # Initialize model via Bouguer approximation
    P = np.zeros(nprisms)
    P = createInitialPrisms(nobs, P, udist, anomObs, prismWidth, prismDensity)

    # Compute initial forward response
    XOBS = anomObs[:, 0]
    YFIT = functn(P, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, udist)

    # Initialize BFGS state
    state = {
        'Hk': np.eye(nprisms),
        'bbeta': np.zeros(nprisms)
    }

    it = 0
    convCont = 0
    chisq1 = None

    while True:
        chisqr, chisqA, tvValue, YFIT, P, flambda, state = QNiterate(
            nprisms, P, anomObs, YFIT,
            xminPrisms, prismWidth, prismDensity,
            dcType, dcBeta, alpha, betaTV,
            flambda, it + 1, udist, state
        )
        it += 1

        print(f"Iter {it:4d} | Obj = {chisqr:.6f} | RMS = {chisqA:.6f} | TV = {tvValue:.6f}")

        # Termination: Marquardt gave up
        if flambda < 0:
            print("Marquardt damping exceeded limit — stopping.")
            break

        # Convergence check
        if it > 1:
            conv = abs(chisq1 - chisqr) / chisq1
            if conv <= usrTOL:
                convCont += 1
                if convCont >= 5:
                    print("Convergence reached.")
                    break
            else:
                convCont = 0

        if it == maxIter:
            print("Maximum number of iterations reached.")
            break

        chisq1 = chisqr

    return P


if __name__ == '__main__':
    import scipy.io
    import os

    # Load observed data from anomObs.mat (MATLAB format)
    # Expected variable 'anomObs': shape (N, 2) — columns: [x, gz]
    mat_path = os.path.join(os.path.dirname(__file__),
                            '..', '..', 'TVGRAV1D', 'anomObs.mat')
    mat = scipy.io.loadmat(mat_path)
    anomObs = mat['anomObs'].astype(float)

    # Parameters matching exemplo1.m from MATLAB
    P = inv_1d(
        anomObs,
        maxIter=400,
        flambda=0.001,
        usrTOL=0.0001,
        alpha=1.0,
        betaTV=0.001,
        nprisms=100,
        xminPrisms=0.0,
        xmaxPrisms=200.0,
        prismDensity=-0.3,
        dcType=0,
        dcBeta=0.0,
        udist=0,
        noise_sigma=1 / 80
    )

    P_result = -1.0 * P  # Sign convention (as in exemplo1.m)

    # Plot result
    import matplotlib.pyplot as plt
    x_prisms = np.linspace(0, 200, len(P_result))
    plt.figure(figsize=(10, 4))
    plt.step(x_prisms, P_result, where='post', color='r', linewidth=2)
    plt.fill_between(x_prisms, P_result, step='post', color='red', alpha=0.3)
    plt.gca().invert_yaxis()
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (km)')
    plt.title('1D TV Gravity Inversion — Final Model')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
