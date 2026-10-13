"""
Calibração da fração f_δ do limiar de Huber nos casos sintéticos.

δ = f_δ × profundidade máxima da estimativa de Bouguer. Para cada f_δ, roda a
curva L de Huber (curva_l.py) nos dois sintéticos, toma o modelo do canto e
calcula o erro RMS de profundidade dividido pela profundidade máxima
verdadeira. A fração escolhida é a de menor erro relativo médio.

Uso (a partir da pasta python/; demora: 5 frações × 2 casos × 18 inversões):

    python experimentos/calibrar_delta.py
    python experimentos/calibrar_delta.py --fracoes 0.025 0.05
"""

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import numpy as np  # noqa: E402

from casos import kael  # noqa: E402
from experimentos import curva_l as cl  # noqa: E402

SINTETICOS = ["graben_isolado", "duas_subbacias"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fracoes", type=float, nargs="+", default=[0.01, 0.025, 0.05, 0.10, 0.20])
    a = p.parse_args()
    erros = {f: [] for f in a.fracoes}
    for nome in SINTETICOS:
        caso = kael.carregar(nome)
        for f in a.fracoes:
            delta = cl.delta_huber(caso, f)
            lambdas = cl.grade_lambdas(cl.centro_grade(nome, "huber", delta))
            res, _, _, i = cl.curva_l(nome, "huber", lambdas, delta)
            e = float(np.sqrt(np.mean((res[i]["m"] - caso.modelo_verdadeiro) ** 2)))
            erros[f].append(e / caso.modelo_verdadeiro.max())
            print(f"{nome}: f = {f:g}  δ = {delta:.4g} km  λ = {lambdas[i]:.4g}  erro RMS = {e:.4f} km")
    media = {f: np.mean(v) for f, v in erros.items()}
    for f, v in media.items():
        print(f"f = {f:g}: erro relativo médio = {100 * v:.2f}%")
    print(f"Fração escolhida: {min(media, key=media.get):g}")


if __name__ == "__main__":
    main()
