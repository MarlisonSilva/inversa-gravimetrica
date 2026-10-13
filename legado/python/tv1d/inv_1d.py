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
    Inversão gravimétrica 1D com regularização por Variação Total (TV).

    Estima a espessura de prismas 2D justapostos (relevo do embasamento)
    minimizando

        T(P) = ||g_obs - F(P)||^2 + alpha * TV(P)

    onde F é a modelagem direta e TV penaliza a variação total do modelo,
    preservando descontinuidades abruptas. A minimização usa Quasi-Newton
    (BFGS) com amortecimento de Marquardt.

    Parâmetros
    ----------
    anomObs : np.ndarray, shape (N, 2)
        Dados observados — colunas: [posição_x, gz].
        gz = -9999.0 sinaliza observação faltante.
    maxIter : int
        Número máximo de iterações.
    flambda : float
        Amortecimento de Marquardt inicial.
    usrTOL : float
        Tolerância de convergência sobre a variação relativa da função objetivo.
    alpha : float
        Peso da regularização TV.
    betaTV : float
        Suavizador que torna a TV diferenciável em torno de zero.
    nprisms : int
        Número de prismas do modelo.
    xminPrisms, xmaxPrisms : float
        Início e fim da discretização ao longo do perfil.
    prismDensity : float
        Contraste de densidade dos prismas (g/cm³).
    dcType : int
        Tipo de contraste (0=constante, 1=hiperbólico — não implementado).
    dcBeta : float
        Fator de decaimento (usado se dcType=1).
    udist : int
        Unidade de distância (0=km, 1=m).
    noise_sigma : float
        Desvio padrão do ruído gaussiano relativo a somar aos dados (0 = sem ruído).

    Retorna
    -------
    P : np.ndarray
        Modelo final — espessura de cada prisma.
    """
    prismWidth = (xmaxPrisms - xminPrisms) / nprisms
    nobs = len(anomObs)

    if noise_sigma > 0.0:
        noise = noise_sigma * np.random.randn(nobs)
        anomObs = anomObs.copy()
        anomObs[:, 1] = anomObs[:, 1] + anomObs[:, 1] * noise

    # Modelo inicial: aproximação de placa de Bouguer a partir dos dados
    P = np.zeros(nprisms)
    P = createInitialPrisms(nobs, P, udist, anomObs, prismWidth, prismDensity)

    # Resposta gravimétrica do modelo inicial
    XOBS = anomObs[:, 0]
    YFIT = functn(P, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, udist)

    # Estado do BFGS entre iterações: Hessiana aproximada (começa em I)
    # e gradiente da iteração anterior
    state = {
        'Hk': np.eye(nprisms),
        'bbeta': np.zeros(nprisms)
    }

    it = 0
    convCont = 0
    chisq1 = None

    while True:
        # Um passo Quasi-Newton: calcula o passo com Hk + lambda, aumenta
        # lambda até a função objetivo cair e atualiza Hk (BFGS)
        chisqr, chisqA, tvValue, YFIT, P, flambda, state = QNiterate(
            nprisms, P, anomObs, YFIT,
            xminPrisms, prismWidth, prismDensity,
            dcType, dcBeta, alpha, betaTV,
            flambda, it + 1, udist, state
        )
        it += 1

        print(f"Iter {it:4d} | Obj = {chisqr:.6f} | RMS = {chisqA:.6f} | TV = {tvValue:.6f}")

        # lambda estourou o limite: nenhum passo reduz mais a função objetivo
        if flambda < 0:
            print("Marquardt damping exceeded limit — stopping.")
            break

        # Convergência: variação relativa pequena por 5 iterações seguidas
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

    # Semente fixa: o ruído (e o modelo invertido) é o mesmo a cada execução
    np.random.seed(7)

    # 1) Modelo verdadeiro sintético: três blocos de profundidades diferentes
    #    com bordas abruptas (caso em que a TV é vantajosa). Profundidades em km.
    nprisms = 100
    xminPrisms = 0.0
    xmaxPrisms = 200.0
    prismWidth = (xmaxPrisms - xminPrisms) / nprisms
    prismDensity = -0.3
    udist = 0  # 0 = km

    model_true = np.zeros(nprisms)
    for i in range(nprisms):
        x = i * prismWidth + prismWidth / 2.0
        if 20 <= x <= 60:
            model_true[i] = 20.0
        elif 80 <= x <= 100:
            model_true[i] = 10.0
        elif 120 <= x <= 180:
            model_true[i] = 30.0

    # 2) Dados observados: modelagem direta do modelo verdadeiro + ruído
    XOBS = np.linspace(0, xmaxPrisms, 51)
    grav_obs = functn(model_true, XOBS, xminPrisms, prismWidth,
                      prismDensity, dcType=0, dcBeta=0.0, unit=udist)

    noise_sigma = 1 / 80  # ruído relativo de 1,25%
    noise = noise_sigma * np.random.randn(len(grav_obs))
    grav_obs_noisy = grav_obs + grav_obs * noise

    anomObs = np.column_stack([XOBS, grav_obs_noisy])

    # 3) Inversão (noise_sigma=0: o ruído já foi somado acima)
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
        noise_sigma=0.0
    )

    # 4) Anomalia calculada a partir do modelo invertido
    grav_calc = functn(P, XOBS, xminPrisms, prismWidth,
                       prismDensity, dcType=0, dcBeta=0.0, unit=udist)

    # 5) Comparação: verdadeiro vs invertido
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
