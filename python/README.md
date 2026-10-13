# Inversão gravimétrica do relevo do embasamento com PSO

Código do artigo: estimativa do relevo do embasamento de uma bacia sedimentar
2D a partir de dados gravimétricos, minimizando

    Φ(m) = Φ_d(m) + λ Φ_r(m),   Φ_d(m) = (1/N) ‖g_obs − g(m)‖₂²

com Otimização por Enxame de Partículas (PSO) e três regularizações:
suavidade (SV), Variação Total (TV) e Huber.

O código antigo (quasi-Newton com TV, MATLAB e a tradução para Python) está
em `../legado/` e não é usado aqui.

## Estrutura

```
python/
├── gravinv/                  pacote principal
│   ├── forward.py            problema direto (prismas 2D, topo na superfície)
│   ├── regularizacao.py      S, SV, TV e Huber
│   ├── objetivo.py           Φ(m) = Φ_d + λ Φ_r
│   ├── pso_serial.py         PSO serial (e o laço comum às duas versões)
│   ├── pso_paralelo.py       PSO com avaliação das partículas em processos
│   └── graficos.py           figuras de ajuste, modelo e convergência
├── casos/
│   ├── kael.py               exemplos de Araújo (2021)
│   └── dados/                dados do graben e da Ponte do Poema
├── experimentos/
│   ├── rodar_exemplo.py      inverte um caso com uma regularização
│   ├── comparar_tempos.py    serial × paralelo (tempo, speedup)
│   ├── curva_l.py            escolha de λ pela curva L (por continuação)
│   ├── calibrar_delta.py     calibração da fração f_δ do limiar de Huber
│   └── parametros.py         λ e f_δ adotados no artigo
├── testes/test_gravinv.py
└── resultados/               saída dos experimentos (criada ao rodar)
```

## Como rodar

A partir desta pasta (`python/`), com o venv ativo:

```
pip install -r requirements.txt
python testes/test_gravinv.py
python experimentos/rodar_exemplo.py --caso graben_isolado --reg tv
python experimentos/rodar_exemplo.py --caso ponte_do_poema --reg huber --processos 2
python experimentos/rodar_exemplo.py --caso graben_isolado --reg sv --repeticoes 10
python experimentos/rodar_exemplo.py --caso duas_subbacias --reg tv --inicializacao uniforme
python experimentos/comparar_tempos.py --caso graben_isolado --reg tv --processos 2 4
python experimentos/curva_l.py --caso graben_isolado --reg huber --fracao-delta 0.025
python experimentos/calibrar_delta.py
```

Casos: `graben_isolado`, `duas_subbacias`, `ponte_do_poema`.
Regularizações: `sv`, `tv`, `huber`.

## Correspondência com o artigo

| Artigo | Código |
|---|---|
| Modelo interpretativo, topos na superfície, p_j = profundidade da base | `forward.ModeloDireto` (`z_topo = 0`) |
| g_i = Σ_j f_i(Δρ₀, p_j, r_ij) | `ModeloDireto.calcular` |
| Φ_d = (1/N)‖g_obs − g(m)‖₂² | `FuncaoObjetivo.residuo` |
| Matriz S, (Sm)_k = p_{k+1} − p_k | `regularizacao.matriz_suavidade`, `variacoes` |
| Φ_SV, Φ_TV, Φ_H (ρ_δ = x² ou 2δ\|x\| − δ²) | `regularizacao.suavidade`, `variacao_total`, `huber` |
| Equações de velocidade e posição, r1 e r2 vetores | `pso_serial.atualizar_enxame` |
| Ciclo do PSO (etapas 1 a 6) | `pso_serial.executar_laco` |
| Paralelização da avaliação das partículas | `pso_paralelo.executar_pso_paralelo` |

## Escolhas registradas no artigo

- Parâmetros do PSO (padrão): K = 300 partículas, 1000 iterações, w = 0,7298,
  c1 = c2 = 1,49618, velocidades iniciais nulas, |v| ≤ 5% do intervalo de busca.
- Posições iniciais: em torno da estimativa da placa de Bouguer (padrão) ou
  uniformes nos limites (seção "Inicialização do enxame").
- Limites: partícula que sai dos limites volta ao limite, com a componente da
  velocidade zerada. Limites por caso em `casos/kael.py`.
- Parada: número máximo de iterações (opcionalmente, estagnação do gbest).
- Paralelização: `multiprocessing` (biblioteca padrão), processos criados por
  "spawn"; o enxame é dividido em blocos iguais, um por processo; o processo
  principal reúne os valores e atualiza o enxame. Serial e paralelo dão
  resultados idênticos com a mesma semente.
- λ: curva L (`curva_l.py`), 18 valores de λ (3 por década), construída por
  continuação (do maior para o menor λ, cada enxame partindo da solução
  anterior); pontos dominados descartados; canto = máximo da curvatura de
  splines de suavização (grau 4, s = 0,03) de log Φ_d e log Φ_r em log λ.
- δ (Huber): δ = f_δ × profundidade máxima da estimativa de Bouguer, com
  f_δ = 2,5% calibrado nos sintéticos (`calibrar_delta.py`).
- Valores adotados: `experimentos/parametros.py`.

## Diferenças em relação ao código legado

1. Forward: o termo da base usava `-h·(θ4 − θ2)` no lugar de `-(z_t + h)·(θ4 − θ2)`.
   O erro era pequeno com topo a 0,5 m, mas cresce com o topo (12% com topo
   a 200 m). Agora o topo fica na superfície, como no artigo, e a fórmula é
   estável com observações sobre as bordas dos prismas.
2. Φ_d agora divide por N (como no artigo e no TCC do Kael).
3. TV é a norma L1 exata; a versão legada era suavizada com β para o quasi-Newton.
4. Derivadas e Jacobiano não existem mais: o PSO não usa derivadas. A
   estimativa de Bouguer (que no legado era o modelo inicial do quasi-Newton)
   virou uma opção de inicialização do enxame.
5. Constante gravitacional: CODATA 2018 (6,67430e-11), no lugar de 6,672e-11.

## Inicialização do enxame

Por padrão, o enxame começa em torno de uma estimativa grosseira do relevo pela placa de Bouguer
(`forward.estimativa_bouguer`): em cada observação, h = g / (2*pi*G*drho),
interpolado para os centros dos prismas. Cada partícula recebe
m0 * (1 + 0,3 * N(0, 1)), limitada aos limites de busca. Com
`--inicializacao uniforme`, as posições são sorteadas dentro dos limites.

Motivo: no caso das duas sub-bacias (80 parâmetros, bacias estreitas e
profundas), a inicialização uniforme deixa o PSO preso em mínimos locais.
O modelo das duas sub-bacias não está no projeto e foi digitalizado da
Figura 6a do TCC do Kael (aproximado).
