"""
Exemplos de Araújo (2021), "Inversão Gravimétrica de Relevo do Embasamento de
Bacia Sedimentar utilizando Algoritmo Genético Paralelizado" (TCC, UFRN).

Três casos, com os parâmetros do TCC:

| Caso             | Seção | Prismas | Extensão | Δρ (g/cm³) | Prof. máx. |
|------------------|-------|---------|----------|------------|------------|
| graben_isolado   | 3.1   | 60      | 60 km    | -0,25      | 1,85 km    |
| duas_subbacias   | 3.2   | 80      | 80 km    | -0,25      | 3,5 km     |
| ponte_do_poema   | 3.3   | 45      | 22,5 m   | -2,3       | ~3,3 m     |

Origem dos dados (detalhes em casos/dados/*.txt):
- graben_isolado: modelo e dados originais do Kael (matlab/modeloGI.mat e
  matlab/gobsGI.mat), com os dados convertidos para mGal.
- duas_subbacias: o modelo NÃO está no projeto. As profundidades abaixo foram
  DIGITALIZADAS à mão da Figura 6a do TCC e são aproximadas. Os dados são
  gerados com o forward deste pacote mais ruído gaussiano (σ = 0,1 mGal,
  como no texto do TCC), com semente fixa.
- ponte_do_poema: dados reais (TVGRAV1D/Dados/Real/Poema). O TCC fala em
  "21,5 m" de extensão, que é o alcance das observações (0 a 21,5 m); a
  discretização segue o arquivo de entrada de Lima (2009): 45 prismas de
  0,5 m, de 0 a 22,5 m.

Os limites de busca do PSO não constam no TCC; os valores abaixo são
escolhas deste projeto (cerca de 1,5x a profundidade máxima esperada).
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from gravinv.forward import ModeloDireto

DADOS = Path(__file__).resolve().parent / "dados"


@dataclass
class Caso:
    nome: str
    modelo_direto: ModeloDireto
    g_obs: np.ndarray
    modelo_verdadeiro: np.ndarray | None   # None para dados reais
    limite_inferior: float
    limite_superior: float
    unidade: str
    descricao: str


def _ler(nome):
    return np.loadtxt(DADOS / nome, comments="#")


def graben_isolado():
    d = _ler("graben_isolado.txt")
    md = ModeloDireto(d[:, 0], 0.0, 60.0, 60, -0.25, unidade="km")
    return Caso("graben_isolado", md, d[:, 1], d[:, 2], 0.0, 3.0, "km",
                "Graben isolado (Araújo, 2021, seção 3.1)")


# Pontos (x em km, profundidade em km) lidos da Figura 6a do TCC.
# Paredes verticais são representadas por dois pontos muito próximos.
_SUBBACIAS_PONTOS = np.array([
    (0.0, 0.00), (2.0, 0.00), (4.0, 0.25), (5.0, 0.32), (13.0, 0.35),
    (15.0, 0.30), (15.01, 3.00), (25.0, 3.00), (25.01, 0.50),
    (28.0, 0.40), (36.0, 0.35), (40.0, 0.25), (44.0, 0.00), (46.0, 0.25),
    (48.0, 0.32), (54.0, 0.35), (56.0, 0.80), (57.5, 2.00), (59.0, 3.00),
    (61.5, 3.50), (63.0, 3.20), (64.0, 3.00), (64.01, 0.50), (66.0, 0.40),
    (69.0, 0.30), (70.0, 0.00), (80.0, 0.00),
])


def modelo_duas_subbacias():
    """Profundidades (km) dos 80 prismas, digitalizadas da Figura 6a (aprox.)."""
    centros = np.arange(80) + 0.5
    return np.round(np.interp(centros, _SUBBACIAS_PONTOS[:, 0],
                              _SUBBACIAS_PONTOS[:, 1]), 3)


def duas_subbacias(sigma_ruido=0.1, semente=8):
    x_obs = np.arange(1.0, 81.0)
    md = ModeloDireto(x_obs, 0.0, 80.0, 80, -0.25, unidade="km")
    verdadeiro = modelo_duas_subbacias()
    rng = np.random.default_rng(semente)
    g_obs = md.calcular(verdadeiro) + sigma_ruido * rng.standard_normal(x_obs.size)
    return Caso("duas_subbacias", md, g_obs, verdadeiro, 0.0, 5.0, "km",
                "Duas sub-bacias (Araújo, 2021, seção 3.2; modelo digitalizado)")


def ponte_do_poema():
    d = _ler("ponte_do_poema.txt")
    md = ModeloDireto(d[:, 0], 0.0, 22.5, 45, -2.3, unidade="m")
    return Caso("ponte_do_poema", md, d[:, 1], None, 0.0, 5.0, "m",
                "Ponte do Poema, UFPA (Araújo, 2021, seção 3.3; dados reais)")


CASOS = {
    "graben_isolado": graben_isolado,
    "duas_subbacias": duas_subbacias,
    "ponte_do_poema": ponte_do_poema,
}


def carregar(nome):
    if nome not in CASOS:
        raise ValueError(f"caso desconhecido: {nome}. Opções: {list(CASOS)}")
    return CASOS[nome]()
