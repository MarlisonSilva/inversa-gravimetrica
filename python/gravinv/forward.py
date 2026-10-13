"""
Problema direto: anomalia gravimétrica de um conjunto de prismas 2D justapostos.

Modelo interpretativo (seção "Modelo Interpretativo" do artigo):
- M prismas retangulares 2D, justapostos, de mesma largura;
- topos de todos os prismas na mesma superfície (z = z_topo, por padrão 0);
- contraste de densidade Δρ constante e conhecido;
- o parâmetro p_j é a profundidade da BASE do j-ésimo prisma.

Convenção: eixo z positivo para baixo; observações em z = 0.

A atração vertical de um prisma de largura [x_e, x_d] e profundidades
[z_t, z_b] num ponto (x_o, 0) é (Telford et al., 1990; equivalente à
eq. A.1 do Apêndice A de Araújo, 2021):

    g = 2 G Δρ [ L(a, z_b) - L(a, z_t) + L(b, z_t) - L(b, z_b)
                 - z_b (θ(b, z_b) - θ(a, z_b))
                 + z_t (θ(b, z_t) - θ(a, z_t)) ]

com a = x_o - x_e, b = x_o - x_d, L(u, z) = (u/2) ln(u² + z²) e
θ(u, z) = arctan(u / z). Escrita assim, a fórmula fica finita mesmo quando
a observação está sobre a borda de um prisma e o topo está na superfície
(z_t = 0): nesses casos os termos u ln|u| e z_t θ tendem a zero.

Correção em relação ao código legado (grav2DPrism.m / grav2DPrism.py):
o termo do arco-tangente da base usava "-h" no lugar de "-(z_t + h)",
o que só é exato quando o topo está em z = 0.
"""

import numpy as np

# Constante gravitacional (CODATA 2018), em m³ kg⁻¹ s⁻²
G = 6.67430e-11
# 1 m/s² = 1e5 mGal
SI_PARA_MGAL = 1e5

_ESCALA = {"km": 1000.0, "m": 1.0}


def _termo_log(u, z):
    """L(u, z) = (u/2) ln(u² + z²), com L = 0 quando u = 0."""
    r2 = u * u + z * z
    seguro = np.where(r2 > 0.0, r2, 1.0)
    return np.where(u != 0.0, 0.5 * u * np.log(seguro), 0.0)


def _termo_angulo(u, z):
    """θ(u, z) = arctan(u / z) para z ≥ 0 (arctan2 evita divisão por zero)."""
    return np.arctan2(u, z)


def estimativa_bouguer(modelo_direto, g_obs):
    """
    Estimativa grosseira do relevo pela aproximação da placa de Bouguer.

    Em cada observação, a anomalia é tratada como se viesse de uma placa
    horizontal infinita de espessura h: g = 2*pi*G*drho*h, logo
    h = g / (2*pi*G*drho). Os valores de h são interpolados para os centros
    dos prismas e limitados a h >= 0.

    Serve como ponto de partida para o enxame do PSO (opção
    `modelo_inicial` de ConfiguracaoPSO). Não é usada como resultado.

    Retorna as profundidades (M,) na unidade do modelo direto.
    """
    md = modelo_direto
    g_obs = np.asarray(g_obs, dtype=float)
    ordem = np.argsort(md.x_obs)
    # 2*pi*G*drho em mGal por metro de espessura
    mgal_por_metro = 2.0 * np.pi * G * md.delta_rho * 1000.0 * SI_PARA_MGAL
    h = g_obs[ordem] / mgal_por_metro / _ESCALA[md.unidade]
    h_prismas = np.interp(md.centros, md.x_obs[ordem], h)
    return np.clip(h_prismas, 0.0, None)


class ModeloDireto:
    """
    Calcula g(m), a anomalia gravimétrica (mGal) do modelo m = [p_1, ..., p_M].

    Parâmetros
    ----------
    x_obs : array (N,)
        Posições das observações ao longo do perfil (na unidade `unidade`).
    x_min, x_max : float
        Início e fim da discretização em prismas.
    n_prismas : int
        Número de prismas M.
    delta_rho : float
        Contraste de densidade bacia − embasamento, em g/cm³ (negativo).
    unidade : {"km", "m"}
        Unidade de todas as distâncias (posições e profundidades).
    z_topo : float
        Profundidade do topo dos prismas. O artigo assume topos na
        superfície (0). Mantido como parâmetro apenas por generalidade.
    """

    def __init__(self, x_obs, x_min, x_max, n_prismas, delta_rho,
                 unidade="km", z_topo=0.0):
        if unidade not in _ESCALA:
            raise ValueError("unidade deve ser 'km' ou 'm'")
        self.x_obs = np.asarray(x_obs, dtype=float)
        self.n_obs = self.x_obs.size
        self.n_prismas = int(n_prismas)
        self.bordas = np.linspace(x_min, x_max, self.n_prismas + 1)
        self.largura = (x_max - x_min) / self.n_prismas
        self.centros = 0.5 * (self.bordas[:-1] + self.bordas[1:])
        self.delta_rho = float(delta_rho)
        self.unidade = unidade
        self.z_topo = float(z_topo)

        # 2 G Δρ (SI) × escala de comprimento × conversão para mGal.
        # Os termos entre colchetes têm dimensão de comprimento.
        self._constante = (2.0 * G * self.delta_rho * 1000.0
                           * _ESCALA[unidade] * SI_PARA_MGAL)

        # Distâncias horizontais observação → bordas esquerda (a) e direita (b),
        # matrizes N × M. Não dependem do modelo: calculadas uma única vez.
        self._a = self.x_obs[:, None] - self.bordas[None, :-1]
        self._b = self.x_obs[:, None] - self.bordas[None, 1:]

        # Parte da fórmula que depende só do topo (fixo): soma sobre os prismas.
        zt = self.z_topo
        topo = (-_termo_log(self._a, zt) + _termo_log(self._b, zt)
                + zt * (_termo_angulo(self._b, zt) - _termo_angulo(self._a, zt)))
        self._parcela_topo = topo.sum(axis=1)  # (N,)

    def calcular(self, m):
        """
        Anomalia calculada g(m) em mGal.

        m : array (M,) com um modelo, ou (K, M) com K modelos (um por linha,
            por exemplo as K partículas do PSO). Retorna (N,) ou (K, N).
        """
        m = np.asarray(m, dtype=float)
        um_modelo = m.ndim == 1
        if um_modelo:
            m = m[None, :]
        if m.shape[1] != self.n_prismas:
            raise ValueError(f"o modelo deve ter {self.n_prismas} profundidades,"
                             f" recebeu {m.shape[1]}")

        zb = m[:, None, :]  # (K, 1, M): profundidade da base de cada prisma
        a, b = self._a[None], self._b[None]
        base = (_termo_log(a, zb) - _termo_log(b, zb)
                - zb * (_termo_angulo(b, zb) - _termo_angulo(a, zb)))
        g = self._constante * (base.sum(axis=2) + self._parcela_topo)
        return g[0] if um_modelo else g
