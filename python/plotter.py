import matplotlib.pyplot as plt

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
