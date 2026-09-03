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
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

    from plotter import plot_inversion_comparison

    # Semente fixa: garante que o ruído — e portanto o modelo invertido —
    # seja o mesmo a cada execução deste exemplo.
    np.random.seed(7)

    # ============================================================
    # 1) Definir o MODELO VERDADEIRO (sintético) para comparação
    # ============================================================
    nprisms = 100
    xminPrisms = 0.0
    xmaxPrisms = 200.0
    prismWidth = (xmaxPrisms - xminPrisms) / nprisms
    prismDensity = -0.3
    udist = 0  # 0 = km

    # Modelo verdadeiro: bacia sedimentar tipo graben
    # Profundidades em km (valores positivos para baixo)
    model_true = np.zeros(nprisms)
    for i in range(nprisms):
        x = i * prismWidth + prismWidth / 2.0
        if 20 <= x <= 60:
            model_true[i] = 20.0   # bloco profundo esquerdo
        elif 80 <= x <= 100:
            model_true[i] = 10.0   # bloco raso central
        elif 120 <= x <= 180:
            model_true[i] = 30.0   # bloco profundo direito

    # ============================================================
    # 2) Gerar dados OBSERVADOS a partir do modelo verdadeiro
    # ============================================================
    XOBS = np.linspace(0, xmaxPrisms, 51)
    grav_obs = functn(model_true, XOBS, xminPrisms, prismWidth,
                      prismDensity, dcType=0, dcBeta=0.0, unit=udist)

    # Adicionar ruído gaussiano (~5%)
    noise_sigma = 1 / 80
    noise = noise_sigma * np.random.randn(len(grav_obs))
    grav_obs_noisy = grav_obs + grav_obs * noise

    # anomObs no formato esperado: (N, 2) — [x, gz]
    anomObs = np.column_stack([XOBS, grav_obs_noisy])

    # ============================================================
    # 3) Executar a INVERSÃO
    # ============================================================
    P = inv_1d(
        anomObs,
        maxIter=400,
        flambda=0.001,
        usrTOL=0.0001,
        alpha=1.0,
        betaTV=0.001,
        nprisms=nprisms,
        xminPrisms=xminPrisms,
        xmaxPrisms=xmaxPrisms,
        prismDensity=prismDensity,
        dcType=0,
        dcBeta=0.0,
        udist=udist,
        noise_sigma=0.0  # ruído já foi adicionado acima
    )

    # ============================================================
    # 4) Calcular a anomalia a partir do modelo invertido
    # ============================================================
    grav_calc = functn(P, XOBS, xminPrisms, prismWidth,
                       prismDensity, dcType=0, dcBeta=0.0, unit=udist)

    # ============================================================
    # 5) Plotar comparação: verdadeiro vs invertido
    # ============================================================
    x_prisms = np.arange(nprisms) * prismWidth

    plot_inversion_comparison(
        x_obs=XOBS,
        grav_obs=grav_obs_noisy,
        grav_calc=grav_calc,
        x_prisms_true=x_prisms,
        model_true=model_true,
        x_prisms_inv=x_prisms,
        model_inv=P,
        title='Inversão Gravimétrica 1D — Variação Total'
    )
