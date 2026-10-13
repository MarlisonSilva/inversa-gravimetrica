"""Figuras: ajuste dos dados, modelo estimado e convergência do PSO."""

import matplotlib.pyplot as plt
import numpy as np


def _degraus(bordas, profundidades):
    """Coordenadas para desenhar os prismas como degraus."""
    x = np.repeat(bordas, 2)[1:-1]
    z = np.repeat(profundidades, 2)
    return x, z


def figura_inversao(caso, modelo_estimado, g_calc, titulo="", modelos_individuais=None):
    """
    Dois painéis: anomalia observada × calculada (topo) e relevo do
    embasamento verdadeiro × estimado (base).
    """
    md = caso.modelo_direto
    u = caso.unidade
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.4]})
    if titulo:
        fig.suptitle(titulo, fontsize=13)

    ax1.plot(md.x_obs, caso.g_obs, "o", ms=4, color="0.2", label="Observada")
    ax1.plot(md.x_obs, g_calc, "-", lw=1.8, color="tab:red", label="Calculada")
    ax1.set_ylabel("Anomalia (mGal)")
    ax1.legend(loc="lower right")
    ax1.grid(alpha=0.3)

    if modelos_individuais is not None:
        for m in modelos_individuais:
            ax2.plot(*_degraus(md.bordas, m), lw=0.6, color="tab:red", alpha=0.25)
    if caso.modelo_verdadeiro is not None:
        ax2.plot(*_degraus(md.bordas, caso.modelo_verdadeiro), lw=2, color="0.2",
                 label="Verdadeiro")
    ax2.plot(*_degraus(md.bordas, modelo_estimado), lw=2, color="tab:red",
             label="Estimado")
    ax2.set_xlabel(f"Distância horizontal ({u})")
    ax2.set_ylabel(f"Profundidade ({u})")
    ax2.set_ylim(bottom=max(caso.limite_superior * 0.8,
                            1.1 * float(np.max(modelo_estimado))), top=0)
    ax2.legend(loc="lower right")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def figura_convergencia(historicos, titulo="Convergência do PSO"):
    """historicos: dict {rótulo: array com a melhor Φ por iteração}."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for rotulo, h in historicos.items():
        ax.plot(np.arange(len(h)), h, lw=1.6, label=rotulo)
    ax.set_yscale("log")
    ax.set_xlabel("Iteração")
    ax.set_ylabel("Melhor Φ(m)")
    ax.set_title(titulo)
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    return fig
