"""
Escolha do parâmetro de regularização λ pela curva L, como em Araújo (2021),
que usa o método da curva L de Hansen (1999).

Para cada λ de uma grade logarítmica, roda-se uma inversão (PSO) e guarda-se
o par (Φ_d, Φ_r) da solução. Em escala log-log, esses pontos formam um "L";
o λ escolhido é o do canto, tomado como o ponto de maior curvatura. A
curvatura é calculada sobre splines de suavização ajustadas a log Φ_d e
log Φ_r em função de log λ (ver curvatura_suavizada).

Continuação em λ: a grade é percorrida do maior para o menor λ. O primeiro
enxame parte da estimativa de Bouguer; cada um dos seguintes parte da solução
do λ anterior, que entra no enxame como uma das partículas. Sem isso, para λ
pequeno o PSO muitas vezes não chega ao mínimo em 1000 iterações e a curva
fica com pontos fora de ordem. Pontos que ainda assim ficam dominados (Φ_d e
Φ_r maiores que os de outro ponto) são descartados antes de achar o canto.

Uso (a partir da pasta python/):

    python experimentos/curva_l.py --caso graben_isolado --reg tv
    python experimentos/curva_l.py --caso duas_subbacias --reg huber --fracao-delta 0.05

Saída em python/resultados/curva_l/: tabela (.txt), figura (.png) e .npz.
"""

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import numpy as np  # noqa: E402

from casos import kael  # noqa: E402
from experimentos import parametros  # noqa: E402
from gravinv.forward import estimativa_bouguer  # noqa: E402
from gravinv.objetivo import FuncaoObjetivo  # noqa: E402
from gravinv.pso_serial import ConfiguracaoPSO, executar_pso  # noqa: E402


def delta_huber(caso, fracao):
    """δ = fração × profundidade máxima da estimativa de Bouguer."""
    m0 = estimativa_bouguer(caso.modelo_direto, caso.g_obs)
    return fracao * float(m0.max())


def centro_grade(nome_caso, reg, delta=None):
    """Centro da grade de λ (valores de parametros.LAMBDA_REFERENCIA)."""
    ref = parametros.LAMBDA_REFERENCIA[nome_caso]
    if reg == "huber":
        # Para δ pequeno, Φ_H ≈ 2δ·TV; centraliza a grade no λ equivalente.
        return ref["tv"] / (2.0 * delta)
    return ref[reg]


def grade_lambdas(centro, k_min=-11, k_max=6):
    """
    λ = centro × 10^(k/3), k = k_min, ..., k_max: três valores por década.
    O padrão (18 valores, de centro × 10^-3,67 a centro × 10^2) cobre a parte
    "vertical" da curva (λ pequeno) e a parte "horizontal" (λ grande).
    """
    return centro * 10.0 ** (np.arange(k_min, k_max + 1) / 3.0)


def inverter(nome_caso, reg, lam, delta=None, semente=0, inicializacao="bouguer",
             n_particulas=300, n_iteracoes=1000, modelo_partida=None):
    """
    Uma inversão; devolve modelo, Φ, Φ_d, Φ_r, histórico e tempo.

    modelo_partida: se dado, o enxame parte dele (e o inclui como partícula),
    em vez da estimativa de Bouguer. Usado na continuação da curva L.
    """
    caso = kael.carregar(nome_caso)
    funcao = FuncaoObjetivo(caso.modelo_direto, caso.g_obs, reg, lam, delta)
    m0, incluir = None, False
    if modelo_partida is not None:
        m0, incluir = np.asarray(modelo_partida, float), True
    elif inicializacao == "bouguer":
        m0 = estimativa_bouguer(caso.modelo_direto, caso.g_obs)
    cfg = ConfiguracaoPSO(caso.limite_inferior, caso.limite_superior,
                          n_particulas=n_particulas, n_iteracoes=n_iteracoes,
                          semente=semente, modelo_inicial=m0,
                          incluir_modelo_inicial=incluir)
    r = executar_pso(funcao, cfg)
    m = r.melhor_modelo
    return {"m": m, "phi": r.melhor_valor, "phi_d": funcao.residuo(m),
            "phi_r": funcao.termo_regularizacao(m), "hist": r.historico,
            "tempo": r.tempo_total}


def curvatura_menger(phi_d, phi_r):
    """
    Curvatura (com sinal) da curva L em log10, por três pontos consecutivos.
    Pontos ordenados por λ crescente. As extremidades recebem NaN.
    Positiva onde a curva dobra no sentido do canto do "L".
    """
    x = np.log10(np.asarray(phi_d, float))
    y = np.log10(np.asarray(phi_r, float))
    k = np.full(x.size, np.nan)
    for i in range(1, x.size - 1):
        ax, ay = x[i] - x[i - 1], y[i] - y[i - 1]
        bx, by = x[i + 1] - x[i], y[i + 1] - y[i]
        cx, cy = x[i + 1] - x[i - 1], y[i + 1] - y[i - 1]
        den = np.hypot(ax, ay) * np.hypot(bx, by) * np.hypot(cx, cy)
        k[i] = 2.0 * (ax * by - ay * bx) / den if den > 0 else np.nan
    return k


def nao_dominados(phi_d, phi_r):
    """
    Máscara dos pontos não dominados da curva L.

    Com o mínimo global de Φ para cada λ, Φ_d cresce e Φ_r decresce com λ.
    Um ponto com Φ_d e Φ_r maiores que os de outro ponto (dominado) indica
    que o PSO não chegou ao mínimo para aquele λ; ele é descartado antes de
    procurar o canto.
    """
    phi_d, phi_r = np.asarray(phi_d, float), np.asarray(phi_r, float)
    mascara = np.ones(phi_d.size, bool)
    for i in range(phi_d.size):
        outros = np.arange(phi_d.size) != i
        if np.any(outros & (phi_d <= phi_d[i]) & (phi_r <= phi_r[i])
                  & ((phi_d < phi_d[i]) | (phi_r < phi_r[i]))):
            mascara[i] = False
    return mascara


def curvatura_suavizada(lambdas, phi_d, phi_r, suavizacao=0.03, n=2000):
    """
    Curvatura da curva L calculada sobre splines de suavização (grau 4)
    ajustadas a log10 Φ_d e log10 Φ_r em função de log10 λ:

        κ(t) = (x' y'' - x'' y') / (x'^2 + y'^2)^(3/2),  t = log10 λ.

    A suavização reduz o efeito das pequenas oscilações dos pontos, que vêm
    do caráter estocástico do PSO. Devolve (t, κ) numa malha fina de t.
    """
    from scipy.interpolate import UnivariateSpline
    t = np.log10(np.asarray(lambdas, float))
    o = np.argsort(t)
    t = t[o]
    x = np.log10(np.asarray(phi_d, float))[o]
    y = np.log10(np.asarray(phi_r, float))[o]
    sx = UnivariateSpline(t, x, k=4, s=suavizacao)
    sy = UnivariateSpline(t, y, k=4, s=suavizacao)
    tt = np.linspace(t[0], t[-1], n)
    x1, x2 = sx.derivative(1)(tt), sx.derivative(2)(tt)
    y1, y2 = sy.derivative(1)(tt), sy.derivative(2)(tt)
    return tt, (x1 * y2 - x2 * y1) / (x1 ** 2 + y1 ** 2) ** 1.5


def indice_canto(lambdas, phi_d, phi_r, suavizacao=0.03):
    """
    Índice do λ do canto da curva L: o ponto da grade (entre os não
    dominados) mais próximo do máximo da curvatura suavizada.
    """
    lambdas = np.asarray(lambdas, float)
    validos = np.flatnonzero(nao_dominados(phi_d, phi_r))
    tt, k = curvatura_suavizada(lambdas[validos], np.asarray(phi_d)[validos],
                                np.asarray(phi_r)[validos], suavizacao)
    t_canto = tt[int(np.argmax(k))]
    return int(validos[np.argmin(np.abs(np.log10(lambdas[validos]) - t_canto))])


def curva_l(nome_caso, reg, lambdas, delta=None, semente=0, ao_terminar=None, **kw):
    """
    Roda a grade de λ por continuação (do maior para o menor λ) e devolve os
    resultados na ordem de `lambdas`, Φ_d, Φ_r e o índice do canto.

    ao_terminar(j, resultado): chamada opcional após cada λ (para salvar).
    """
    lambdas = np.asarray(lambdas, float)
    res = [None] * lambdas.size
    partida = None
    for j in np.argsort(lambdas)[::-1]:
        res[j] = inverter(nome_caso, reg, lambdas[j], delta, semente,
                          modelo_partida=partida, **kw)
        partida = res[j]["m"]
        if ao_terminar is not None:
            ao_terminar(j, res[j])
    phi_d = np.array([r["phi_d"] for r in res])
    phi_r = np.array([r["phi_r"] for r in res])
    return res, phi_d, phi_r, indice_canto(lambdas, phi_d, phi_r)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--caso", choices=list(kael.CASOS), default="graben_isolado")
    p.add_argument("--reg", choices=["sv", "tv", "huber"], default="tv")
    p.add_argument("--fracao-delta", type=float, default=parametros.FRACAO_DELTA,
                   help="δ de Huber como fração da prof. máxima de Bouguer")
    p.add_argument("--semente", type=int, default=0)
    a = p.parse_args()

    caso = kael.carregar(a.caso)
    delta = delta_huber(caso, a.fracao_delta) if a.reg == "huber" else None
    lambdas = grade_lambdas(centro_grade(a.caso, a.reg, delta))
    res, phi_d, phi_r, i = curva_l(a.caso, a.reg, lambdas, delta, a.semente)

    saida = RAIZ / "resultados" / "curva_l"
    saida.mkdir(parents=True, exist_ok=True)
    base = saida / (f"{a.caso}_{a.reg}" + (f"_d{a.fracao_delta:g}" if delta else ""))
    linhas = [f"# {caso.descricao} | {a.reg}" + (f" | δ = {delta:.4g} {caso.unidade}" if delta else ""),
              "# lambda\tphi_d\tphi_r"]
    linhas += [f"{lam:.6g}\t{d:.6g}\t{r:.6g}" for lam, d, r in zip(lambdas, phi_d, phi_r)]
    linhas.append(f"# canto: lambda = {lambdas[i]:.6g}")
    Path(f"{base}.txt").write_text("\n".join(linhas) + "\n", encoding="utf-8")
    np.savez(f"{base}.npz", lambdas=lambdas, phi_d=phi_d, phi_r=phi_r, canto=i,
             modelos=np.array([r["m"] for r in res]), delta=delta if delta else np.nan)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 4))
    ok = nao_dominados(phi_d, phi_r)
    ax.loglog(phi_d[ok], phi_r[ok], "o-")
    ax.loglog(phi_d[~ok], phi_r[~ok], "x", color="0.5")
    ax.loglog(phi_d[i], phi_r[i], "o", ms=12, mfc="none", mec="r")
    ax.set_xlabel("Φ_d")
    ax.set_ylabel("Φ_r")
    ax.set_title(f"Curva L: {a.caso}, {a.reg}, λ = {lambdas[i]:.3g}")
    fig.savefig(f"{base}.png", dpi=150, bbox_inches="tight")
    print("\n".join(linhas))
    print(f"\nArquivos salvos em {saida}")


if __name__ == "__main__":
    main()
