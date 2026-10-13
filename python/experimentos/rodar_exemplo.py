"""
Inverte um exemplo do Kael com PSO e uma regularização.

Uso (a partir da pasta python/):

    python experimentos/rodar_exemplo.py --caso graben_isolado --reg tv
    python experimentos/rodar_exemplo.py --caso ponte_do_poema --reg huber --processos 2
    python experimentos/rodar_exemplo.py --caso graben_isolado --reg sv --repeticoes 10
    python experimentos/rodar_exemplo.py --caso duas_subbacias --reg tv --inicializacao uniforme

Com --repeticoes N, o PSO roda N vezes com sementes diferentes e o modelo
final é a média dos N modelos, como no TCC do Kael (média de 10 inversões).
--processos 0 (padrão) usa o PSO serial; --processos K usa o paralelo.
Por padrão, o enxame começa em torno da estimativa da placa de Bouguer
(--inicializacao bouguer); --inicializacao uniforme sorteia as posições
dentro dos limites. λ e δ padrão: experimentos/parametros.py (curva L).

Saída em python/resultados/: figura da inversão (.png), curva de
convergência (.png) e os modelos (.npz).
"""

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import matplotlib  # noqa: E402

import numpy as np  # noqa: E402

from casos import kael  # noqa: E402
from experimentos import parametros  # noqa: E402
from experimentos.curva_l import delta_huber  # noqa: E402
from gravinv.forward import estimativa_bouguer  # noqa: E402
from gravinv.objetivo import FuncaoObjetivo  # noqa: E402
from gravinv.pso_serial import ConfiguracaoPSO, executar_pso  # noqa: E402
from gravinv.pso_paralelo import executar_pso_paralelo  # noqa: E402

NOMES_REG = {"sv": "Suavidade (SV)", "tv": "Variação Total (TV)", "huber": "Huber"}


def argumentos():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--caso", choices=list(kael.CASOS), default="graben_isolado")
    p.add_argument("--reg", choices=["sv", "tv", "huber"], default="tv")
    p.add_argument("--lam", type=float, help="λ (padrão: experimentos/parametros.py)")
    p.add_argument("--delta", type=float,
                   help="δ de Huber (padrão: FRACAO_DELTA × prof. máx. de Bouguer)")
    p.add_argument("--particulas", type=int, default=300)
    p.add_argument("--iteracoes", type=int, default=1000)
    p.add_argument("--repeticoes", type=int, default=1)
    p.add_argument("--semente", type=int, default=0)
    p.add_argument("--processos", type=int, default=0, help="0 = PSO serial")
    p.add_argument("--inicializacao", choices=["uniforme", "bouguer"], default="bouguer")
    p.add_argument("--mostrar", action="store_true", help="abre as figuras na tela")
    return p.parse_args()


def main():
    a = argumentos()
    if not a.mostrar:
        matplotlib.use("Agg")
    from gravinv.graficos import figura_convergencia, figura_inversao

    caso = kael.carregar(a.caso)
    lam = a.lam if a.lam is not None else parametros.LAMBDA[a.caso][a.reg]
    delta = None
    if a.reg == "huber":
        delta = a.delta if a.delta is not None else delta_huber(caso, parametros.FRACAO_DELTA)
    funcao = FuncaoObjetivo(caso.modelo_direto, caso.g_obs, a.reg, lam, delta)

    print(f"{caso.descricao}")
    print(f"Regularização: {NOMES_REG[a.reg]}  λ = {lam:g}"
          + (f"  δ = {delta:g} {caso.unidade}" if delta else ""))
    modo = "serial" if a.processos == 0 else f"paralelo ({a.processos} processos)"
    print(f"PSO {modo}: {a.particulas} partículas, {a.iteracoes} iterações,"
          f" {a.repeticoes} repetição(ões), inicialização {a.inicializacao}\n")
    m0 = None
    if a.inicializacao == "bouguer":
        m0 = estimativa_bouguer(caso.modelo_direto, caso.g_obs)

    modelos, historicos, tempos = [], {}, []
    for r in range(a.repeticoes):
        cfg = ConfiguracaoPSO(caso.limite_inferior, caso.limite_superior,
                              n_particulas=a.particulas, n_iteracoes=a.iteracoes,
                              semente=a.semente + r, modelo_inicial=m0)
        if a.processos == 0:
            res = executar_pso(funcao, cfg)
        else:
            res = executar_pso_paralelo(funcao, cfg, n_processos=a.processos)
        modelos.append(res.melhor_modelo)
        historicos[f"semente {cfg.semente}"] = res.historico
        tempos.append(res.tempo_total)
        print(f"  repetição {r + 1}: Φ = {res.melhor_valor:.5g}"
              f"  ({res.n_iteracoes} iterações, {res.tempo_total:.1f} s)")

    modelos = np.array(modelos)
    medio = modelos.mean(axis=0)
    g_calc = caso.modelo_direto.calcular(medio)
    phi_d = funcao.residuo(medio)
    phi_r = funcao.termo_regularizacao(medio)
    print(f"\nModelo final ({'média das repetições' if a.repeticoes > 1 else 'melhor partícula'}):")
    print(f"  Φ_d = {phi_d:.5g} mGal²   (RMS = {np.sqrt(phi_d):.4g} mGal)")
    print(f"  Φ_r = {phi_r:.5g}   Φ = {phi_d + lam * phi_r:.5g}")
    print(f"  profundidade máxima = {medio.max():.3f} {caso.unidade}")
    if caso.modelo_verdadeiro is not None:
        erro = np.sqrt(np.mean((medio - caso.modelo_verdadeiro) ** 2))
        print(f"  RMS do erro de profundidade = {erro:.3f} {caso.unidade}")
    print(f"  tempo médio por inversão = {np.mean(tempos):.1f} s")

    saida = RAIZ / "resultados"
    saida.mkdir(exist_ok=True)
    base = saida / f"{a.caso}_{a.reg}_{a.inicializacao}"
    titulo = f"{caso.descricao}\n{NOMES_REG[a.reg]}, λ = {lam:g}, inicialização {a.inicializacao}"
    figura_inversao(caso, medio, g_calc, titulo,
                    modelos if a.repeticoes > 1 else None).savefig(f"{base}.png", dpi=150)
    figura_convergencia(historicos).savefig(f"{base}_convergencia.png", dpi=150)
    np.savez(f"{base}.npz", modelos=modelos, modelo_medio=medio, g_calc=g_calc,
             lam=lam, delta=delta if delta else np.nan, tempos=np.array(tempos))
    print(f"\nArquivos salvos em {saida}")
    if a.mostrar:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
