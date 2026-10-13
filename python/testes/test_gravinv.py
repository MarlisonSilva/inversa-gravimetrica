"""
Testes do pacote gravinv.

Rodar a partir da pasta python/:
    python testes/test_gravinv.py        (sem dependências extras)
    python -m pytest testes              (se tiver pytest)
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import numpy as np  # noqa: E402

from gravinv import regularizacao as reg  # noqa: E402
from gravinv.forward import G, ModeloDireto, estimativa_bouguer  # noqa: E402
from gravinv.objetivo import FuncaoObjetivo  # noqa: E402
from gravinv.pso_serial import ConfiguracaoPSO, executar_pso, inicializar_enxame  # noqa: E402
from gravinv.pso_paralelo import executar_pso_paralelo  # noqa: E402


def _integral_numerica(drho, x_e, x_d, z_t, z_b, x_o, n=2000):
    """g de um prisma por soma de linhas de massa (km, g/cm³ → mGal)."""
    xs = np.linspace(x_e, x_d, n + 1)
    zs = np.linspace(z_t, z_b, n + 1)
    X, Z = np.meshgrid(0.5 * (xs[:-1] + xs[1:]) - x_o, 0.5 * (zs[:-1] + zs[1:]))
    soma = (Z / (X ** 2 + Z ** 2)).sum() * (x_d - x_e) / n * (z_b - z_t) / n
    return 2 * G * drho * 1000 * soma * 1000 * 1e5


def test_forward_contra_integracao_numerica():
    for z_t, z_b, x_o in [(0.2, 1.2, 0.3), (0.5, 1.0, 2.0), (0.0, 1.0, 3.0)]:
        md = ModeloDireto([x_o], 0.0, 1.0, 1, -0.25, "km", z_topo=z_t)
        g = md.calcular([z_b])[0]
        ref = _integral_numerica(-0.25, 0.0, 1.0, z_t, z_b, x_o)
        assert abs(g - ref) < 1e-3 * abs(ref), (z_t, z_b, x_o, g, ref)


def test_forward_placa_de_bouguer():
    # prisma muito largo ≈ placa infinita: 2πGΔρh
    md = ModeloDireto([0.0], -5000.0, 5000.0, 1, -0.25, "km")
    esperado = 2 * np.pi * G * (-250.0) * 1000.0 * 1e5
    assert abs(md.calcular([1.0])[0] - esperado) < 1e-3 * abs(esperado)


def test_forward_topo_na_superficie_e_bordas_sem_nan():
    # observações exatamente sobre as bordas dos prismas
    md = ModeloDireto(np.arange(0.0, 11.0), 0.0, 10.0, 10, -0.25, "km")
    g = md.calcular(np.linspace(0.0, 2.0, 10))
    assert np.all(np.isfinite(g))
    # prisma de espessura nula não gera anomalia
    assert np.allclose(md.calcular(np.zeros(10)), 0.0)


def test_forward_unidades_km_e_m():
    km = ModeloDireto([0.3, 1.5], 0.0, 1.0, 2, -2.3, "km")
    m = ModeloDireto([300.0, 1500.0], 0.0, 1000.0, 2, -2.3, "m")
    assert np.allclose(km.calcular([0.4, 0.8]), m.calcular([400.0, 800.0]))


def test_forward_varios_modelos_de_uma_vez():
    md = ModeloDireto(np.linspace(0, 20, 15), 0.0, 20.0, 8, -0.25, "km")
    M = np.random.default_rng(0).random((5, 8)) * 2
    em_lote = md.calcular(M)
    um_a_um = np.array([md.calcular(m) for m in M])
    assert np.allclose(em_lote, um_a_um)


def test_regularizacoes():
    m = np.array([0.0, 0.0, 1.0, 1.05, 3.0])
    d = m @ reg.matriz_suavidade(5).T
    assert np.allclose(d, reg.variacoes(m))
    assert np.isclose(reg.suavidade(m), np.sum(d ** 2))
    assert np.isclose(reg.variacao_total(m), np.sum(np.abs(d)))
    # Huber: x² até δ, 2δ|x| − δ² acima; contínua em |x| = δ
    delta = 0.5
    esperado = sum(x * x if abs(x) <= delta else 2 * delta * abs(x) - delta ** 2 for x in d)
    assert np.isclose(reg.huber(m, delta), esperado)
    eps = 1e-9
    assert np.isclose(reg.funcao_huber(delta - eps, delta), reg.funcao_huber(delta + eps, delta))
    # δ grande → igual à SV; δ pequeno → ≈ 2δ·TV
    assert np.isclose(reg.huber(m, 100.0), reg.suavidade(m))
    assert np.isclose(reg.huber(m, 1e-6), 2e-6 * reg.variacao_total(m), rtol=1e-3)


def test_residuo_divide_por_n():
    md = ModeloDireto([1.0, 2.0, 3.0], 0.0, 4.0, 4, -0.25, "km")
    m = np.array([0.5, 1.0, 1.0, 0.5])
    g = md.calcular(m)
    f = FuncaoObjetivo(md, g + np.array([1.0, -1.0, 2.0]), "tv", 0.0)
    assert np.isclose(f.residuo(m), (1 + 1 + 4) / 3)


def test_pso_encontra_o_minimo_de_uma_quadratica():
    class Quadratica:
        n_parametros = 5
        def __call__(self, X):
            return np.sum((np.atleast_2d(X) - 1.5) ** 2, axis=-1)
    cfg = ConfiguracaoPSO(-5.0, 5.0, n_particulas=40, n_iteracoes=300, semente=1)
    r = executar_pso(Quadratica(), cfg)
    assert np.allclose(r.melhor_modelo, 1.5, atol=1e-3)
    assert np.all(np.diff(r.historico) <= 0)  # o gbest nunca piora


def test_pso_paralelo_igual_ao_serial():
    md = ModeloDireto(np.arange(1.0, 21.0), 0.0, 20.0, 20, -0.25, "km")
    verdadeiro = np.where((np.arange(20) > 5) & (np.arange(20) < 14), 1.5, 0.3)
    f = FuncaoObjetivo(md, md.calcular(verdadeiro), "huber", 0.01, delta=0.1)
    cfg = ConfiguracaoPSO(0.0, 3.0, n_particulas=30, n_iteracoes=50, semente=7)
    s = executar_pso(f, cfg)
    p = executar_pso_paralelo(f, cfg, n_processos=2)
    assert np.array_equal(s.melhor_modelo, p.melhor_modelo)
    assert np.array_equal(s.historico, p.historico)


def test_estimativa_bouguer_recupera_placa():
    # placa larga e uniforme: a estimativa deve devolver a espessura da placa
    md = ModeloDireto(np.arange(-50.0, 51.0), -5000.0, 5000.0, 1000, -0.25, "km")
    m = np.full(1000, 1.2)
    h = estimativa_bouguer(md, md.calcular(m))
    centrais = np.abs(md.centros) < 40
    assert np.allclose(h[centrais], 1.2, rtol=1e-3)
    # mesma resposta em metros
    md_m = ModeloDireto(np.arange(-50.0, 51.0) * 1000, -5e6, 5e6, 1000, -0.25, "m")
    h_m = estimativa_bouguer(md_m, md_m.calcular(m * 1000))
    assert np.allclose(h_m[centrais], 1200.0, rtol=1e-3)


def test_inicializacao_com_modelo_inicial():
    md = ModeloDireto(np.arange(1.0, 21.0), 0.0, 20.0, 20, -0.25, "km")
    verdadeiro = np.where((np.arange(20) > 5) & (np.arange(20) < 14), 1.5, 0.3)
    g = md.calcular(verdadeiro)
    m0 = estimativa_bouguer(md, g)
    cfg = ConfiguracaoPSO(0.0, 3.0, n_particulas=30, n_iteracoes=40, semente=3,
                          modelo_inicial=m0)
    grande = ConfiguracaoPSO(0.0, 3.0, n_particulas=400, modelo_inicial=m0)
    x, _ = inicializar_enxame(grande, 20, np.random.default_rng(0))
    assert np.all((x >= 0.0) & (x <= 3.0))
    # posições espalhadas em torno de m0 (perturbação relativa de 30%)
    assert np.allclose(np.median(x, axis=0), m0, rtol=0.1)
    f = FuncaoObjetivo(md, g, "tv", 0.01)
    s = executar_pso(f, cfg)
    p = executar_pso_paralelo(f, cfg, n_processos=2)
    assert np.array_equal(s.melhor_modelo, p.melhor_modelo)


def test_incluir_modelo_inicial_como_particula():
    m0 = np.linspace(0.1, 2.0, 20)
    cfg = ConfiguracaoPSO(0.0, 3.0, n_particulas=10, modelo_inicial=m0,
                          incluir_modelo_inicial=True)
    x, _ = inicializar_enxame(cfg, 20, np.random.default_rng(0))
    assert np.array_equal(x[0], m0)
    assert not np.array_equal(x[1], m0)


def test_curva_l_pontos_dominados_e_canto():
    from experimentos.curva_l import indice_canto, nao_dominados
    # ponto "ruim" (índice 2): Φ_d e Φ_r maiores que os do índice 3
    phi_d = np.array([0.010, 0.011, 0.05, 0.013, 0.05, 0.3, 2.0])
    phi_r = np.array([100.0, 40.0, 30.0, 10.0, 8.0, 6.0, 5.0])
    ok = nao_dominados(phi_d, phi_r)
    assert list(ok) == [True, True, False, True, True, True, True]
    # curva L sintética: ramo vertical (Φ_d ≈ 0,01) e horizontal (Φ_r ≈ 1),
    # canto em λ = 1
    lam = 10.0 ** np.arange(-4, 4.01, 1 / 3)
    phi_d = 0.01 * (1 + lam ** 2)
    phi_r = 1 + 1 / lam ** 2
    assert abs(np.log10(lam[indice_canto(lam, phi_d, phi_r)])) <= 1 / 3 + 1e-9


def test_pso_respeita_limites():
    md = ModeloDireto(np.arange(1.0, 11.0), 0.0, 10.0, 10, -0.25, "km")
    f = FuncaoObjetivo(md, md.calcular(np.full(10, 2.9)), "sv", 0.0)
    r = executar_pso(f, ConfiguracaoPSO(0.0, 1.0, n_particulas=20, n_iteracoes=30, semente=0))
    assert np.all((r.melhor_modelo >= 0.0) & (r.melhor_modelo <= 1.0))


if __name__ == "__main__":
    testes = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in testes:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(testes)} testes passaram.")
