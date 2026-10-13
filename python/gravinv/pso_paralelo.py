"""
Otimização por Enxame de Partículas (PSO), versão paralela.

Segue a seção "Paralelização do PSO" do artigo: só a avaliação das
partículas é paralela. Em cada iteração, o enxame é dividido em blocos
(um por processo); cada processo calcula o problema direto e a função
objetivo das partículas do seu bloco; os valores voltam ao processo
principal, que atualiza pbest, gbest, velocidades e posições de forma
sequencial.

A inicialização, a atualização do enxame e o critério de parada são as
mesmas funções de pso_serial.py. Com a mesma semente, o resultado é igual
ao da versão serial; muda apenas o tempo de execução.

Recurso: multiprocessing (biblioteca padrão do Python), com processos
criados por "spawn", que é o único modo disponível no Windows. Por isso,
scripts que chamam executar_pso_paralelo precisam do bloco
`if __name__ == "__main__":`.
"""

import multiprocessing as mp
import os
import time

import numpy as np

from .pso_serial import executar_laco

# Função objetivo de cada processo trabalhador. É enviada uma única vez,
# quando o processo é criado, e não a cada iteração.
_FUNCAO = None


def _inicializar_processo(funcao):
    global _FUNCAO
    _FUNCAO = funcao


def _avaliar_bloco(bloco):
    return _FUNCAO(bloco)


def executar_pso_paralelo(funcao, cfg, n_processos=None, n_parametros=None):
    """
    PSO com avaliação das partículas distribuída entre processos.

    funcao : objeto chamável Φ(X) que aceita (K, M) e devolve (K,).
             Precisa ser "picklable" (FuncaoObjetivo é).
    n_processos : número de processos trabalhadores (padrão: núcleos
                  lógicos da máquina, limitado ao número de partículas).
    """
    if n_parametros is None:
        n_parametros = funcao.n_parametros
    if n_processos is None:
        n_processos = os.cpu_count() or 1
    n_processos = max(1, min(int(n_processos), cfg.n_particulas))

    inicio = time.perf_counter()
    contexto = mp.get_context("spawn")
    with contexto.Pool(n_processos, initializer=_inicializar_processo,
                       initargs=(funcao,)) as pool:
        tempo_criar_processos = time.perf_counter() - inicio

        def avaliar(x):
            blocos = np.array_split(x, n_processos)
            return np.concatenate(pool.map(_avaliar_bloco, blocos))

        resultado = executar_laco(avaliar, n_parametros, cfg)

    resultado.tempo_total = time.perf_counter() - inicio
    resultado.extras["n_processos"] = n_processos
    resultado.extras["tempo_criar_processos"] = tempo_criar_processos
    return resultado
