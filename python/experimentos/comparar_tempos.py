"""
Compara o tempo do PSO serial com o do PSO paralelo.

Para cada semente, roda o serial e o paralelo com a mesma configuração,
confere que os dois chegaram ao mesmo modelo (o algoritmo é o mesmo; só
muda a forma de avaliar as partículas) e mede o tempo. Ao final mostra o
tempo médio, o speedup (serial / paralelo) e a redução percentual.

Uso (a partir da pasta python/):

    python experimentos/comparar_tempos.py --caso graben_isolado --reg tv --processos 2 4
"""

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import numpy as np  # noqa: E402

from casos import kael  # noqa: E402
from experimentos import parametros  # noqa: E402
from experimentos.curva_l import delta_huber  # noqa: E402
from gravinv.forward import estimativa_bouguer  # noqa: E402
from gravinv.objetivo import FuncaoObjetivo  # noqa: E402
from gravinv.pso_serial import ConfiguracaoPSO, executar_pso  # noqa: E402
from gravinv.pso_paralelo import executar_pso_paralelo  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--caso", choices=list(kael.CASOS), default="graben_isolado")
    p.add_argument("--reg", choices=["sv", "tv", "huber"], default="tv")
    p.add_argument("--processos", type=int, nargs="+", default=[2])
    p.add_argument("--particulas", type=int, default=300)
    p.add_argument("--iteracoes", type=int, default=1000)
    p.add_argument("--repeticoes", type=int, default=3)
    p.add_argument("--inicializacao", choices=["uniforme", "bouguer"], default="bouguer")
    a = p.parse_args()

    caso = kael.carregar(a.caso)
    lam = parametros.LAMBDA[a.caso][a.reg]
    delta = delta_huber(caso, parametros.FRACAO_DELTA) if a.reg == "huber" else None
    funcao = FuncaoObjetivo(caso.modelo_direto, caso.g_obs, a.reg, lam, delta)
    m0 = (estimativa_bouguer(caso.modelo_direto, caso.g_obs)
          if a.inicializacao == "bouguer" else None)

    tempos = {"serial": []} | {f"{n} processos": [] for n in a.processos}
    for s in range(a.repeticoes):
        cfg = ConfiguracaoPSO(caso.limite_inferior, caso.limite_superior,
                              n_particulas=a.particulas, n_iteracoes=a.iteracoes,
                              semente=s, modelo_inicial=m0)
        serial = executar_pso(funcao, cfg)
        tempos["serial"].append(serial.tempo_total)
        linha = f"semente {s}: serial {serial.tempo_total:.1f} s"
        for n in a.processos:
            par = executar_pso_paralelo(funcao, cfg, n_processos=n)
            if not np.array_equal(par.melhor_modelo, serial.melhor_modelo):
                raise RuntimeError("serial e paralelo deram modelos diferentes")
            tempos[f"{n} processos"].append(par.tempo_total)
            linha += f" | {n} proc. {par.tempo_total:.1f} s"
        print(linha, flush=True)

    t_serial = np.mean(tempos["serial"])
    print(f"\n{caso.descricao}, {a.reg.upper()}: {a.particulas} partículas,"
          f" {a.iteracoes} iterações, {a.repeticoes} repetições")
    print(f"{'Execução':<14}{'tempo médio (s)':>16}{'speedup':>10}{'redução':>10}")
    for nome, t in tempos.items():
        tm = np.mean(t)
        print(f"{nome:<14}{tm:>16.1f}{t_serial / tm:>10.2f}{100 * (1 - tm / t_serial):>9.1f}%")
    print("\nResultados idênticos entre serial e paralelo em todas as sementes.")


if __name__ == "__main__":
    main()
