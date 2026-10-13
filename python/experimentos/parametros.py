"""
Parâmetros de regularização dos exemplos do Kael.

LAMBDA_REFERENCIA: valores de uma varredura rápida inicial. Servem só como
centro da grade de λ da curva L (experimentos/curva_l.py).

LAMBDA e FRACAO_DELTA: valores finais, escolhidos pela curva L
(experimentos/curva_l.py) e pela calibração de δ nos sintéticos
(experimentos/calibrar_delta.py). São os usados no artigo.

δ (Huber) = FRACAO_DELTA × profundidade máxima da estimativa de Bouguer do
caso (km nos sintéticos, m na Ponte do Poema).
"""

LAMBDA_REFERENCIA = {
    "graben_isolado": {"sv": 0.1, "tv": 0.02},
    "duas_subbacias": {"sv": 0.03, "tv": 0.015},
    "ponte_do_poema": {"sv": 1e-2, "tv": 1e-2},
}

FRACAO_DELTA = 0.025  # calibrada nos sintéticos (experimentos/calibrar_delta.py)

# Canto da curva L (curva_l.py, grade padrão, semente 0). Para Huber, com
# δ = FRACAO_DELTA × profundidade máxima de Bouguer.
LAMBDA = {
    "graben_isolado": {"sv": 0.1, "tv": 0.09283, "huber": 0.4853},
    "duas_subbacias": {"sv": 6.463e-05, "tv": 0.0003232, "huber": 0.005644},
    "ponte_do_poema": {"sv": 1e-05, "tv": 2.154e-06, "huber": 0.0001766},
}
