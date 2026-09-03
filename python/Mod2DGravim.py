from modGrav1dTOTAL import modGrav1dTOTAL

def Mod2DGravim(m, unit=0):
    """
    Main script/function for 2D gravity modeling.

    Distances are in kilometers, so the anomaly is returned in mGal.
    """
    nprisms = 60  # Number of prisms in the interpretation model.
    xminPrisms = 0  # Start of discretization along profile.
    xmaxPrisms = 60  # End of the discretization along profile.
    prismWidth = 1  # Prisms horizontal width, the same for all prisms.
    prismDensity = -0.25  # g/cm3! Density contrast of the prisms, the same for all prisms.
    nobs = 60  # Mudar se necessario
    xobs = list(range(1, 61))  # Equivalent to 1:60 in MATLAB (values 1 to 60)

    Dgrav = modGrav1dTOTAL(m, xobs, xminPrisms, prismWidth, prismDensity, unit)

    return Dgrav

if __name__ == '__main__':
    from plotter import plot_model_and_anomaly
    
    # Criando um modelo sintético de teste (60 prismas)
    m_teste = [0] * 10 + [5] * 10 + [10] * 20 + [5] * 10 + [0] * 10
    
    # Executa a função principal
    anomalia_calculada = Mod2DGravim(m_teste)
    
    # Define o eixo X (1 a 60 para observações e 0 a 59 para os prismas)
    x_obs = list(range(1, 61))
    x_prisms = list(range(len(m_teste)))
    
    # Plota os resultados
    plot_model_and_anomaly(x_obs, anomalia_calculada, x_prisms, m_teste)
