"""
Termos de regularização Φ_r(m) da seção "Métodos de Regularização" do artigo.

Os três atuam sobre as variações entre profundidades de prismas vizinhos,
(Sm)_k = p_{k+1} − p_k, k = 1, ..., M−1, em que S é a matriz de diferenças
finitas de primeira ordem. Diferem apenas na forma de penalizar (Sm)_k:

    SV    : Σ (Sm)_k²                      (norma L2 ao quadrado)
    TV    : Σ |(Sm)_k|                     (norma L1, sem suavização)
    Huber : Σ ρ_δ((Sm)_k),  ρ_δ(x) = x²            se |x| ≤ δ
                                     2δ|x| − δ²    se |x| > δ

Todas as funções aceitam um modelo (M,) ou vários modelos (K, M) e devolvem
um escalar ou um vetor (K,).

Diferença em relação ao código legado: a TV legada era suavizada,
Σ √((Sm)_k² + β²), porque o quasi-Newton exige derivadas. O PSO não usa
derivadas, então aqui a TV é a norma L1 exata, como no artigo.
"""

import numpy as np


def matriz_suavidade(n_prismas):
    """Matriz S, (M−1) × M, de diferenças finitas de primeira ordem."""
    S = np.zeros((n_prismas - 1, n_prismas))
    k = np.arange(n_prismas - 1)
    S[k, k] = -1.0
    S[k, k + 1] = 1.0
    return S


def variacoes(m):
    """Sm: diferenças entre profundidades consecutivas (equivale a S @ m)."""
    return np.diff(np.asarray(m, dtype=float), axis=-1)


def suavidade(m):
    """Φ_SV(m) = ‖Sm‖₂²."""
    d = variacoes(m)
    return np.sum(d * d, axis=-1)


def variacao_total(m):
    """Φ_TV(m) = ‖Sm‖₁."""
    return np.sum(np.abs(variacoes(m)), axis=-1)


def funcao_huber(x, delta):
    """ρ_δ(x): quadrática para |x| ≤ δ e linear para |x| > δ."""
    if delta <= 0:
        raise ValueError("delta deve ser positivo")
    a = np.abs(x)
    return np.where(a <= delta, a * a, 2.0 * delta * a - delta * delta)


def huber(m, delta):
    """Φ_H(m) = Σ ρ_δ((Sm)_k)."""
    return np.sum(funcao_huber(variacoes(m), delta), axis=-1)
