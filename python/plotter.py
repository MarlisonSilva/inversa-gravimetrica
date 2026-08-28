import matplotlib.pyplot as plt
import numpy as np


def plot_model_and_anomaly(x_obs, Dgrav, x_prisms, m):
    """
    Plota a anomalia gravimétrica e o modelo geológico (em formato de blocos/stairs).
    
    Parâmetros:
    x_obs: Lista de coordenadas X das observações
    Dgrav: Anomalia gravimétrica calculada
    x_prisms: Lista de coordenadas X dos prismas
    m: Modelo geológico (profundidade/espessura dos prismas)
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Gráfico 1: Anomalia Gravimétrica
    ax1.plot(x_obs, Dgrav, 'b.-', label='Anomalia Gravimétrica')
    ax1.set_ylabel('Gravidade (mGal)')
    ax1.set_title('Modelagem Direta Gravimétrica')
    ax1.grid(True)
    ax1.legend()

    # Gráfico 2: Modelo Geológico / Prismas (Gráfico tipo "stairs")
    ax2.step(x_prisms, m, where='post', color='r', linewidth=2, label='Modelo (Prismas)')
    ax2.fill_between(x_prisms, m, step='post', color='red', alpha=0.3)
    ax2.set_xlabel('Distância ao longo do perfil')
    ax2.set_ylabel('Profundidade / Espessura')
    ax2.invert_yaxis() # Profundidade crescendo para baixo
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.show()


def plot_inversion_comparison(x_obs, grav_obs, grav_calc,
                              x_prisms_true, model_true,
                              x_prisms_inv, model_inv,
                              title='Inversão Gravimétrica com Variação Total'):
    """
    Plota a comparação entre dados observados/calculados e modelo verdadeiro/invertido
    em uma única figura com dois painéis.

    Painel superior: Anomalia observada vs calculada.
    Painel inferior: Modelo verdadeiro vs modelo invertido (ambos em formato stairs).

    Parâmetros
    ----------
    x_obs : array-like
        Posições das observações gravimétricas.
    grav_obs : array-like
        Anomalia gravimétrica observada.
    grav_calc : array-like
        Anomalia gravimétrica calculada pelo modelo invertido.
    x_prisms_true : array-like
        Eixo X do modelo verdadeiro (posições dos prismas).
    model_true : array-like
        Espessuras/profundidades do modelo verdadeiro.
    x_prisms_inv : array-like
        Eixo X do modelo invertido.
    model_inv : array-like
        Espessuras/profundidades do modelo invertido.
    title : str
        Título geral da figura.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                                    gridspec_kw={'height_ratios': [1, 1.5]})
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # --- Painel superior: Anomalia Gravimétrica ---
    ax1.plot(x_obs, grav_obs, 'ko', markersize=4, label='Dado Observado')
    ax1.plot(x_obs, grav_calc, 'r-', linewidth=1.5, label='Dado Calculado (Inversão)')
    ax1.set_ylabel('Anomalia Gravimétrica (mGal)')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)

    # --- Painel inferior: Modelo Geológico ---
    # Modelo verdadeiro
    ax2.step(x_prisms_true, model_true, where='post', color='blue',
             linewidth=2, label='Modelo Verdadeiro')
    ax2.fill_between(x_prisms_true, model_true, step='post',
                     color='blue', alpha=0.15)
    # Modelo invertido
    ax2.step(x_prisms_inv, model_inv, where='post', color='red',
             linewidth=2, linestyle='--', label='Modelo Invertido (TV)')
    ax2.fill_between(x_prisms_inv, model_inv, step='post',
                     color='red', alpha=0.15)

    ax2.set_xlabel('Distância ao longo do perfil (km)')
    ax2.set_ylabel('Profundidade (km)')
    ax2.invert_yaxis()
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
