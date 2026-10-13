# Código legado

Material anterior ao código do artigo, guardado só como referência. Nada aqui
é importado pelo código novo em `../python/`.

| Pasta | Conteúdo |
|---|---|
| `TVGRAV1D/` | Inversão quasi-Newton com Variação Total (Lima, 2009), em MATLAB, com os dados dos casos em `Dados/` |
| `matlab/` | Modelagem direta e testes do Kael (inclui `modeloGI.mat` e `gobsGI.mat`, usados em `../python/casos/dados/graben_isolado.txt`) |
| `python/` | Tradução para Python do `TVGRAV1D` (pacote `tv1d`) e da modelagem direta do `matlab/`, com `ARQUITETURA.md` |
| `RELATORIO_CODIGO_PYTHON.pdf` | Relatório sobre a tradução para Python |

Erros conhecidos, mantidos aqui como estavam (corrigidos no código novo):

1. `grav2DPrism` (todas as versões): o termo da base usa `-h*(theta4-theta2)`
   no lugar de `-(z0+h)*(theta4-theta2)`.
2. `dgrav2DPrism` (MATLAB e Python): a derivada está deslocada de meio prisma.
   No `.m`, a linha `x = xobs - x0 + w/2.0;` ficou partida e o `+ w/2.0;`
   solto não tem efeito; a tradução repetiu o erro.

Os scripts de `python/` usam caminhos relativos à própria pasta e continuam
rodando a partir de `legado/python/`.
