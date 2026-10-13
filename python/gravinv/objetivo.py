"""
Função objetivo do problema inverso (seção "Problema Inverso" do artigo):

    Φ(m) = Φ_d(m) + λ Φ_r(m),    Φ_d(m) = (1/N) ‖g_obs − g(m)‖₂²

em que Φ_r é um dos termos de regularização (SV, TV ou Huber).

Diferença em relação ao código legado: o resíduo legado não dividia por N.
Isso muda a escala de λ: valores de λ do código legado não valem aqui.
"""

import numpy as np

from . import regularizacao

REGULARIZACOES = ("sv", "tv", "huber")


class FuncaoObjetivo:
    """
    Φ(m) para um modelo direto, dados observados e regularização escolhida.

    Parâmetros
    ----------
    modelo_direto : gravinv.forward.ModeloDireto
    g_obs : array (N,)
        Anomalias observadas (mGal), nas posições modelo_direto.x_obs.
    regularizacao : {"sv", "tv", "huber"}
    lam : float
        Parâmetro de regularização λ (≥ 0).
    delta : float, opcional
        Limiar δ da função de Huber (mesma unidade das profundidades).
        Obrigatório quando regularizacao = "huber".

    Um objeto desta classe pode ser enviado a outros processos (é
    "picklable"), o que o PSO paralelo usa.
    """

    def __init__(self, modelo_direto, g_obs, regularizacao="tv", lam=0.0,
                 delta=None):
        regularizacao = regularizacao.lower()
        if regularizacao not in REGULARIZACOES:
            raise ValueError(f"regularizacao deve ser uma de {REGULARIZACOES}")
        if regularizacao == "huber" and (delta is None or delta <= 0):
            raise ValueError("a regularização de Huber exige delta > 0")
        if lam < 0:
            raise ValueError("lam deve ser ≥ 0")
        self.modelo_direto = modelo_direto
        self.g_obs = np.asarray(g_obs, dtype=float)
        if self.g_obs.shape != (modelo_direto.n_obs,):
            raise ValueError("g_obs deve ter uma anomalia por observação")
        self.regularizacao = regularizacao
        self.lam = float(lam)
        self.delta = delta
        self.n_parametros = modelo_direto.n_prismas

    def residuo(self, m):
        """Φ_d(m) = (1/N) ‖g_obs − g(m)‖₂² (escalar ou vetor (K,))."""
        r = self.g_obs - self.modelo_direto.calcular(m)
        return np.mean(r * r, axis=-1)

    def termo_regularizacao(self, m):
        """Φ_r(m) para a regularização escolhida."""
        if self.regularizacao == "sv":
            return regularizacao.suavidade(m)
        if self.regularizacao == "tv":
            return regularizacao.variacao_total(m)
        return regularizacao.huber(m, self.delta)

    def __call__(self, m):
        """Φ(m) = Φ_d(m) + λ Φ_r(m). Aceita (M,) ou (K, M)."""
        return self.residuo(m) + self.lam * self.termo_regularizacao(m)
