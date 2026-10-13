"""
Otimização por Enxame de Partículas (PSO), versão serial.

Segue a seção "Otimização por Enxame de Partículas (PSO)" do artigo:

    v_k ← w v_k + c1 r1 ⊙ (x_k^pbest − x_k) + c2 r2 ⊙ (x^gbest − x_k)
    x_k ← x_k + v_k

com r1 e r2 vetores de números aleatórios uniformes em [0, 1] sorteados a
cada iteração para cada partícula (⊙ = produto elemento a elemento).

Ciclo:
    1. gerar K partículas com posições aleatórias dentro dos limites;
    2–3. avaliar a função objetivo de cada partícula (problema direto incluso);
    4. atualizar as melhores posições individuais e a melhor global;
    5. atualizar velocidades e posições;
    6. repetir 2–5 até o critério de parada.

Escolhas de implementação (precisam constar no artigo, nos [PREENCHER]):
- posições iniciais uniformes nos limites (padrão) ou em torno de um modelo
  inicial, como a estimativa da placa de Bouguer (opção modelo_inicial);
- velocidades iniciais nulas;
- velocidade limitada a ±fracao_vmax × (limite_superior − limite_inferior);
- partícula que sai dos limites volta para o limite e tem a componente da
  velocidade daquela direção zerada;
- parada pelo número máximo de iterações ou, se `tolerancia` for dada,
  quando a melhor função objetivo melhorar menos que `tolerancia`
  (relativo) ao longo de `paciencia` iterações.

A versão paralela (pso_paralelo.py) usa exatamente as mesmas funções de
inicialização e atualização deste módulo; só muda a forma de avaliar as
partículas. Com a mesma semente, as duas devolvem o mesmo resultado.
"""

import time
from dataclasses import dataclass, field

import numpy as np


@dataclass
class ConfiguracaoPSO:
    """Parâmetros do PSO."""
    limite_inferior: float | np.ndarray
    limite_superior: float | np.ndarray
    n_particulas: int = 300
    n_iteracoes: int = 1000
    w: float = 0.7298          # peso de inércia (constante)
    c1: float = 1.49618        # coeficiente cognitivo
    c2: float = 1.49618        # coeficiente social
    fracao_vmax: float = 0.05  # |v| ≤ 5% do intervalo de busca
    tolerancia: float | None = None
    paciencia: int = 100
    semente: int | None = None
    # Inicialização: None = posições uniformes nos limites (padrão).
    # Com um modelo m0 (M,), as partículas começam em
    # m0 * (1 + perturbacao * N(0, 1)), limitadas aos limites de busca
    # (por exemplo, m0 = forward.estimativa_bouguer(...)).
    modelo_inicial: np.ndarray | None = None
    perturbacao: float = 0.3
    # Se True, a partícula 0 começa exatamente em modelo_inicial (sem
    # perturbação). Usado na curva L por continuação: o enxame de um λ parte
    # da solução do λ anterior e nunca termina pior que ela.
    incluir_modelo_inicial: bool = False


@dataclass
class ResultadoPSO:
    """Saída do PSO."""
    melhor_modelo: np.ndarray
    melhor_valor: float
    historico: np.ndarray            # melhor Φ ao fim de cada iteração
    n_iteracoes: int
    criterio_parada: str
    tempo_total: float               # segundos
    tempo_avaliacao: float           # segundos gastos avaliando partículas
    extras: dict = field(default_factory=dict)


def limites(cfg, n_parametros):
    """Limites inferior e superior como vetores (M,)."""
    li = np.broadcast_to(np.asarray(cfg.limite_inferior, float), (n_parametros,))
    ls = np.broadcast_to(np.asarray(cfg.limite_superior, float), (n_parametros,))
    if np.any(ls <= li):
        raise ValueError("limite_superior deve ser maior que limite_inferior")
    return li.copy(), ls.copy()


def inicializar_enxame(cfg, n_parametros, rng):
    """
    Etapa 1: posições iniciais e velocidades nulas.

    Sem modelo_inicial: posições uniformes dentro dos limites.
    Com modelo_inicial m0: posições m0 * (1 + perturbacao * N(0, 1)),
    limitadas aos limites de busca.
    """
    li, ls = limites(cfg, n_parametros)
    forma = (cfg.n_particulas, n_parametros)
    if cfg.modelo_inicial is None:
        x = li + rng.random(forma) * (ls - li)
    else:
        m0 = np.asarray(cfg.modelo_inicial, dtype=float)
        if m0.shape != (n_parametros,):
            raise ValueError("modelo_inicial deve ter um valor por parâmetro")
        x = np.clip(m0 * (1.0 + cfg.perturbacao * rng.standard_normal(forma)), li, ls)
        if cfg.incluir_modelo_inicial:
            x[0] = np.clip(m0, li, ls)
    v = np.zeros_like(x)
    return x, v


def atualizar_enxame(x, v, pbest_x, gbest_x, cfg, rng):
    """Etapa 5: novas velocidades e posições (equações do artigo)."""
    li, ls = limites(cfg, x.shape[1])
    vmax = cfg.fracao_vmax * (ls - li)
    r1 = rng.random(x.shape)
    r2 = rng.random(x.shape)
    v = cfg.w * v + cfg.c1 * r1 * (pbest_x - x) + cfg.c2 * r2 * (gbest_x - x)
    v = np.clip(v, -vmax, vmax)
    x = x + v
    fora = (x < li) | (x > ls)
    x = np.clip(x, li, ls)
    v = np.where(fora, 0.0, v)
    return x, v


def atualizar_melhores(x, f, pbest_x, pbest_f):
    """Etapa 4: atualiza pbest (no lugar) e devolve o índice do gbest."""
    melhorou = f < pbest_f
    pbest_x[melhorou] = x[melhorou]
    pbest_f[melhorou] = f[melhorou]
    return int(np.argmin(pbest_f))


def executar_laco(avaliar, n_parametros, cfg):
    """
    Laço do PSO, comum às versões serial e paralela.

    avaliar : função que recebe as posições (K, M) e devolve Φ (K,).
    """
    rng = np.random.default_rng(cfg.semente)
    inicio = time.perf_counter()
    tempo_avaliacao = 0.0

    x, v = inicializar_enxame(cfg, n_parametros, rng)
    t = time.perf_counter()
    f = np.asarray(avaliar(x), dtype=float)
    tempo_avaliacao += time.perf_counter() - t

    pbest_x, pbest_f = x.copy(), f.copy()
    g = int(np.argmin(pbest_f))
    historico = [pbest_f[g]]
    criterio = "número máximo de iterações"

    for it in range(1, cfg.n_iteracoes + 1):
        x, v = atualizar_enxame(x, v, pbest_x, pbest_x[g], cfg, rng)
        t = time.perf_counter()
        f = np.asarray(avaliar(x), dtype=float)
        tempo_avaliacao += time.perf_counter() - t
        g = atualizar_melhores(x, f, pbest_x, pbest_f)
        historico.append(pbest_f[g])

        if cfg.tolerancia is not None and it >= cfg.paciencia:
            antigo = historico[-1 - cfg.paciencia]
            if antigo - historico[-1] <= cfg.tolerancia * abs(antigo):
                criterio = (f"melhora relativa ≤ {cfg.tolerancia:g} em "
                            f"{cfg.paciencia} iterações")
                break

    return ResultadoPSO(
        melhor_modelo=pbest_x[g].copy(),
        melhor_valor=float(pbest_f[g]),
        historico=np.array(historico),
        n_iteracoes=len(historico) - 1,
        criterio_parada=criterio,
        tempo_total=time.perf_counter() - inicio,
        tempo_avaliacao=tempo_avaliacao,
    )


def executar_pso(funcao, cfg, n_parametros=None):
    """
    PSO serial: todas as partículas são avaliadas no processo principal.

    funcao : objeto chamável Φ(X) que aceita (K, M) e devolve (K,),
             por exemplo gravinv.objetivo.FuncaoObjetivo.
    """
    if n_parametros is None:
        n_parametros = funcao.n_parametros
    return executar_laco(funcao, n_parametros, cfg)
