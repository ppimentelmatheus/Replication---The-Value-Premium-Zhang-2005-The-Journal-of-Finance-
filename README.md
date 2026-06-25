# Replicacao em Python de Zhang (2005), "The Value Premium"

Este projeto porta para Python a replicacao computacional do artigo:

> Lu Zhang (2005), "The Value Premium", The Journal of Finance, 60(1), 67-103.

O foco aqui e reproduzir a parte de simulacao do modelo: calibracao, solucao do problema dinamico da firma, equilibrio aproximado a la Krusell-Smith, simulacao de paineis de firmas, carteiras value/growth, figuras e tabelas geradas pelo modelo. A parte empirica com dados historicos ainda nao e o objetivo desta base.

Os arquivos originais em MATLAB/Fortran ficam em `vpCode/`. O port em Python fica no pacote `zhang2005/`.

## Estado Atual

Ja implementado:

- Calibracao mensal do benchmark, portada de `vpCode/CalibCC.m`.
- Rouwenhorst para os processos de produtividade agregada e idiossincratica.
- Taxa livre de risco fechada, portada de `getRfcc.m`.
- Iteracao da funcao valor, portada de `vfi3fcnIEccB.f90`.
- Simulacao de equilibrio, distribuicao estacionaria e painel, a partir dos MEX Fortran originais.
- Construcao de SMB/HML e value premium, portada de `ValPrem.m`.
- Figuras basicas do artigo:
  - Figura 1: custo de ajuste assimetrico.
  - Figura 2: lucratividade value/growth no estilo Fama-French (1995).
  - Figura 3: investimento e custo de ajuste em bons e maus tempos.
  - Figura 4: versao simulada/proxy de HML e value spread contra produtividade agregada.
  - Figura 5: momentos do pricing kernel.
  - Figura B.1: qualidade da agregacao aproximada.
- Tabelas model-only:
  - Parametros do modelo.
  - Momentos agregados simulados.
  - Carteiras ordenadas por book-to-market.
  - Regressoes preditivas internas ao modelo.
- Notebooks gerados a partir dos arquivos Python.

Ainda nao implementado completamente:

- Benchmark completo em escala do artigo com `N = 5000`, horizontes longos e `nkp = 5000`.
- Aceleracao robusta com Numba/F2PY para rodadas grandes.
- Figura 4 estrutural exata com retorno esperado por estado. A versao atual e uma proxy simulada usando HML realizado e value spread.
- Table IV de estatica comparativa completa, que exige multiplas rodadas do modelo com parametros alternativos.
- Tabelas empiricas com dados historicos.

## Estrutura do Projeto

```text
.
|-- vpCode/                         # Codigo MATLAB/Fortran original
|-- zhang2005/                      # Pacote Python da replicacao
|   |-- calibration.py              # Calibracao, Rouwenhorst, taxa livre de risco
|   |-- vfi.py                      # Iteracao da funcao valor
|   |-- simulation.py               # Simulacoes de equilibrio e painel
|   |-- equilibrium.py              # Loop Krusell-Smith
|   |-- portfolios.py               # SMB, HML e value premium
|   |-- ff95.py                     # Port da analise FF95 de lucratividade
|   |-- figure_data.py              # Dados intermediarios para figuras
|   |-- model_tables.py             # Tabelas model-only
|   |-- interpolation.py            # Interpolacao 4D das politicas
|   |-- linalg.py                   # OLS
|   `-- io.py                       # Leitura/escrita MAT/NPZ
|-- scripts/
|   |-- run_calibration.py          # Checagem da calibracao benchmark
|   |-- run_mini_replication.py     # Pipeline reduzido end-to-end
|   |-- run_figure_data.py          # Gera dados para figuras
|   |-- plot_article_figures.py     # Plota figuras em PNG/PDF
|   |-- run_model_tables.py         # Gera tabelas model-only
|   `-- run_value_premium.py        # Calcula ValPrem a partir de um .mat
|-- notebooks/                      # Notebooks gerados a partir dos .py
|-- docs/                           # Notas tecnicas do port
|-- tests/                          # Testes smoke
|-- outputs/                        # Resultados gerados por scripts
|-- requirements.txt
`-- pyproject.toml
```

## Ambiente

Use Python 3.10+.

Crie um ambiente virtual local:

```bash
cd "/home/matheus/mpp/1. FGV EPGE Doutorado/2 ano/6 trimestre/Finanças/Replication"
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Dependencias principais:

- `numpy`
- `scipy`
- `numba`
- `matplotlib`

Observacao: o `.venv/` e ignorado pelo `.gitignore`.

## Validacao Basica

Rode os testes:

```bash
.venv/bin/python -m unittest discover -s tests
```

Resultado esperado:

```text
Ran 5 tests ...
OK
```

Cheque a calibracao:

```bash
.venv/bin/python scripts/run_calibration.py
```

Saida esperada aproximada:

```text
xbar: -5.7037210796
capital grid points: 25
k min/max: 0.0100 / 9.9947
kp grid points: 5000
h grid: [2.75  2.875 3.    3.125 3.25 ]
annualized average Sharpe ratio: 0.405524
```

A calibracao Python foi comparada com `vpCode/Params.mat` nos testes smoke.

## Pipeline Original MATLAB

O pipeline original do Zhang e:

1. `CalibCC.m`
   - calibra parametros;
   - gera grids;
   - salva `Params.mat`.

2. `mainCC.m`
   - resolve a funcao valor;
   - estima a lei de movimento aproximada do preco agregado;
   - salva `vfi3Mat.mat` e `coefIS.mat`.

3. `ssIE.m`
   - simula ate a distribuicao estacionaria de firmas;
   - salva `distrSS.mat`.

4. `panIE.m` ou `panIE900.m`
   - simula paineis de firmas;
   - construi carteiras;
   - calcula momentos, tabelas e figuras.

O pipeline Python atual segue a mesma ordem, mas com uma versao reduzida para depuracao e desenvolvimento.

## Pipeline Python: Replicacao Reduzida

O script principal e:

```bash
.venv/bin/python scripts/run_mini_replication.py \
  --periods 181 \
  --ks-iterations 5 \
  --vfi-iterations 120 \
  --output-dir outputs/replication_for_figures
```

Ele executa:

```text
calibracao reduzida
-> solucao da VFI
-> loop Krusell-Smith
-> simulacao de painel
-> SMB/HML
-> value premium
-> salvamento dos resultados
```

Arquivos gerados:

```text
outputs/replication_for_figures/mini_replication_results.npz
outputs/replication_for_figures/mini_replication_summary.json
```

O `.npz` contem objetos como:

- `Pf`: valor de mercado ex-dividendo das firmas.
- `Bf`: book value/capital das firmas.
- `Df`: dividendos.
- `Rf`: retornos das firmas.
- `Rm`: retorno agregado/value-weighted.
- `srf`: taxa livre de risco simulada.
- `sx`: choque agregado simulado.
- `optK`: politica otima de capital.
- `V`: funcao valor.
- `SMB`, `HML`: fatores simulados.
- `table`: tabela de value premium no formato de `ValPrem.m`.

O `.json` contem diagnosticos legiveis da rodada, como coeficientes da lei agregada, `R2`, erros finais da VFI e momentos agregados.

## Pipeline Python: Replicacao Completa

O script para o benchmark completo e:

```text
scripts/run_full_replication.py
```

Ele separa as fases que no MATLAB aparecem em `mainCC.m`, `ssIE.m` e `panIE.m`:

```text
equilibrio Krusell-Smith
-> distribuicao estacionaria
-> simulacoes repetidas de painel
-> value premium por painel
-> resumo agregado
```

Sempre comece com:

```bash
.venv/bin/python scripts/run_full_replication.py \
  --dry-run \
  --output-dir outputs/full_replication_dry_run
```

Para testar o simulador de painel completo sem resolver a VFI em Python, use a solucao MATLAB ja incluida em `vpCode/`:

```bash
.venv/bin/python scripts/run_full_replication.py \
  --use-matlab-equilibrium \
  --stationary-mat vpCode/distrSS.mat \
  --panel-periods 421 \
  --panel-simulations 20 \
  --save-panel-arrays first \
  --output-dir outputs/full_replication_matlab_equilibrium
```

Para rodar tudo do zero em Python, o comando benchmark seria:

```bash
.venv/bin/python scripts/run_full_replication.py \
  --n-firms 5000 \
  --equilibrium-periods 11000 \
  --stationary-periods 10000 \
  --panel-periods 421 \
  --panel-simulations 20 \
  --nkp 5000 \
  --cutoff 1000 \
  --ks-max-iterations 100 \
  --vfi-max-iterations 10000 \
  --output-dir outputs/full_replication_python
```

Aviso: a versao Python pura ainda nao esta otimizada para esse comando completo. A rota pratica no momento e usar `--use-matlab-equilibrium` para validar as fases de simulacao e, em paralelo, acelerar VFI/simulacao com Numba ou F2PY.

## Figuras do Artigo

Depois de rodar a simulacao base, gere os dados das figuras:

```bash
.venv/bin/python scripts/run_figure_data.py \
  --replication-npz outputs/replication_for_figures/mini_replication_results.npz \
  --output-dir outputs/figures/data
```

Depois plote:

```bash
.venv/bin/python scripts/plot_article_figures.py \
  --data-dir outputs/figures/data \
  --png-dir outputs/figures/png \
  --pdf-dir outputs/figures/pdf
```

Para usar uma rodada completa gerada por `scripts/run_full_replication.py`, use:

```bash
.venv/bin/python scripts/run_figure_data.py \
  --full-output-dir outputs/full_replication \
  --recompute-b1 \
  --output-dir outputs/figures/full_data

.venv/bin/python scripts/plot_article_figures.py \
  --data-dir outputs/figures/full_data \
  --png-dir outputs/figures/full_png \
  --pdf-dir outputs/figures/full_pdf
```

No modo full, figuras que exigem painéis completos usam o primeiro arquivo salvo em
`outputs/full_replication/panels/panel_*_arrays.npz`. Se a rodada foi feita com
`--save-panel-arrays none`, apenas as figuras que dependem de calibracao ou fatores
compactos poderao ser geradas. A Figura B.1 exige
`outputs/full_replication/equilibrium_diagnostics.npz`, salvo pelas rodadas full
feitas com a versao atual do script. Se o arquivo nao existir, use `--recompute-b1`
para reconstruir o diagnostico a partir de `equilibrium_solution.npz`.

Arquivos esperados:

```text
outputs/figures/png/figure1_adjustment_cost.png
outputs/figures/png/figure2_ff95_profitability.png
outputs/figures/png/figure3_investment_scatter.png
outputs/figures/png/figure4_simulated_spreads.png
outputs/figures/png/figure5_pricing_kernel.png
outputs/figures/png/figureB1_aggregation_quality.png

outputs/figures/pdf/figure1_adjustment_cost.pdf
outputs/figures/pdf/figure2_ff95_profitability.pdf
outputs/figures/pdf/figure3_investment_scatter.pdf
outputs/figures/pdf/figure4_simulated_spreads.pdf
outputs/figures/pdf/figure5_pricing_kernel.pdf
outputs/figures/pdf/figureB1_aggregation_quality.pdf
```

Notas sobre as figuras:

- Figura 1 e Figura 5 dependem basicamente da calibracao.
- Figura 2 usa `ff95.py`, port de `FF95.m`.
- Figura 3 usa `Pf`, `Bf`, `In` e `sx` do painel simulado.
- Figura B.1 usa os objetos do loop Krusell-Smith.
- Figura 4 atual e uma proxy simulada, nao ainda a figura estrutural exata do artigo.

## Tabelas Model-Only

Depois de rodar a simulacao base, gere as tabelas:

```bash
.venv/bin/python scripts/run_model_tables.py \
  --replication-npz outputs/replication_for_figures/mini_replication_results.npz \
  --summary-json outputs/replication_for_figures/mini_replication_summary.json \
  --output-dir outputs/tables/model
```

Para gerar tabelas usando a rodada completa:

```bash
.venv/bin/python scripts/run_model_tables.py \
  --full-output-dir outputs/full_replication \
  --output-dir outputs/tables/full_model
```

No modo full, o script agrega os arquivos `panel_*_factors.npz` para a tabela
media de value premium e ratios agregados. Tabelas que exigem `Pf`, `Bf` e `Rf`
usam o primeiro painel completo salvo em `panel_*_arrays.npz`, a menos que a
rodada tenha sido feita com `--save-panel-arrays all`.

A Tabela 2 em `table2_model_moments` usa as definicoes do paper: Sharpe anual,
taxa real anual, retorno value-weighted da industria, volatilidade da industria,
volatilidade media individual, book-to-market industrial e taxas de
investimento/desinvestimento. Momentos auxiliares ficam em
`table2_additional_model_moments`.

Arquivos gerados em Markdown:

```text
outputs/tables/model/markdown/table1_model_parameters.md
outputs/tables/model/markdown/table2_model_moments.md
outputs/tables/model/markdown/table3_bm_portfolios_model.md
outputs/tables/model/markdown/table5_6_predictive_regressions_model.md
```

Arquivos equivalentes em CSV:

```text
outputs/tables/model/csv/table1_model_parameters.csv
outputs/tables/model/csv/table2_model_moments.csv
outputs/tables/model/csv/table3_bm_portfolios_model.csv
outputs/tables/model/csv/table5_6_predictive_regressions_model.csv
```

As tabelas atuais sao apenas do modelo. Elas nao misturam dados empiricos.

## Ordem Recomendada Para Reproduzir Tudo

Use esta sequencia:

```bash
cd "/home/matheus/mpp/1. FGV EPGE Doutorado/2 ano/6 trimestre/Finanças/Replication"

.venv/bin/python -m unittest discover -s tests

.venv/bin/python scripts/run_calibration.py

.venv/bin/python scripts/run_mini_replication.py \
  --periods 181 \
  --ks-iterations 5 \
  --vfi-iterations 120 \
  --output-dir outputs/replication_for_figures

.venv/bin/python scripts/run_figure_data.py \
  --replication-npz outputs/replication_for_figures/mini_replication_results.npz \
  --output-dir outputs/figures/data

.venv/bin/python scripts/plot_article_figures.py \
  --data-dir outputs/figures/data \
  --png-dir outputs/figures/png \
  --pdf-dir outputs/figures/pdf

.venv/bin/python scripts/run_model_tables.py \
  --replication-npz outputs/replication_for_figures/mini_replication_results.npz \
  --summary-json outputs/replication_for_figures/mini_replication_summary.json \
  --output-dir outputs/tables/model
```

## Notebooks

A pasta `notebooks/` contem notebooks gerados a partir dos scripts e modulos Python.

Exemplos:

```text
notebooks/scripts/run_calibration.ipynb
notebooks/scripts/run_value_premium.ipynb
notebooks/zhang2005/calibration.ipynb
notebooks/zhang2005/vfi.ipynb
notebooks/zhang2005/simulation.ipynb
notebooks/zhang2005/equilibrium.ipynb
notebooks/zhang2005/portfolios.ipynb
```

Esses notebooks sao conversoes do codigo, nao notebooks narrativos finais. Para uso em texto de replicacao, o ideal e criar depois notebooks organizados assim:

```text
01_calibracao.ipynb
02_equilibrio_vfi.ipynb
03_simulacao_paineis.ipynb
04_figuras.ipynb
05_tabelas.ipynb
```

## Scripts

### `scripts/run_calibration.py`

Imprime objetos principais da calibracao benchmark.

```bash
.venv/bin/python scripts/run_calibration.py
```

### `scripts/run_mini_replication.py`

Executa o pipeline reduzido end-to-end.

Argumentos uteis:

- `--n-firms`: numero de firmas simuladas.
- `--periods`: numero de meses simulados.
- `--nkp`: numero de candidatos para capital futuro.
- `--nx`, `--nz`: tamanhos dos grids de choques.
- `--ks-iterations`: iteracoes Krusell-Smith.
- `--vfi-iterations`: iteracoes da funcao valor.
- `--initial-coefficients`: arquivo `.mat` com `alp1`-`alp4`; padrao `vpCode/coefIS.mat`.
- `--output-dir`: pasta de saida.

Exemplo maior, ainda experimental:

```bash
.venv/bin/python scripts/run_mini_replication.py \
  --n-firms 500 \
  --periods 421 \
  --ks-iterations 10 \
  --vfi-iterations 250 \
  --output-dir outputs/replication_medium
```

### `scripts/run_figure_data.py`

Gera os `.npz` com dados intermediarios das figuras.

### `scripts/plot_article_figures.py`

Le os `.npz` das figuras e salva PNG/PDF.

### `scripts/run_model_tables.py`

Gera tabelas model-only em CSV e Markdown.

### `scripts/run_value_premium.py`

Calcula a tabela de value premium a partir de um `.mat` que contenha `Pf`, `Rf`, `Bf`, `Rm` e `srf` ou `rf`.

```bash
.venv/bin/python scripts/run_value_premium.py caminho/arquivo.mat
```

## Mapeamento MATLAB -> Python

| MATLAB/Fortran original | Python |
| --- | --- |
| `CalibCC.m` | `zhang2005/calibration.py` |
| `rouwTrans.m` | `zhang2005/calibration.py::rouwenhorst` |
| `CspSimu.m` | `zhang2005/calibration.py::simulate_ar1` |
| `getRfcc.m` | `zhang2005/calibration.py::get_rf_cc` |
| `vfi3fcnIEccB.f90` | `zhang2005/vfi.py` |
| `simIEfcn3.f90` | `zhang2005/simulation.py::simulate_equilibrium` |
| `ssIEfcn.f90` | `zhang2005/simulation.py::simulate_stationary_distribution` |
| `panIEfcn.f90` | `zhang2005/simulation.py::simulate_panel` |
| `mainCC.m` | `zhang2005/equilibrium.py` |
| `ValPrem.m` | `zhang2005/portfolios.py` |
| `FF95.m` | `zhang2005/ff95.py` |
| `lambda_m.m` | `zhang2005/figure_data.py::pricing_kernel_moments` |
| `ols.m` | `zhang2005/linalg.py` |

## Observacoes Computacionais

O benchmark completo do artigo e grande:

```text
N = 5000 firmas
Ts = 11000 ou 12000 meses no equilibrio
nkp = 5000 candidatos de capital futuro
nx = 11
nz = 15
nh = 5
```

A versao Python atual prioriza fidelidade e legibilidade. Para rodar o benchmark completo, sera necessario acelerar principalmente:

- `zhang2005/vfi.py`
- `zhang2005/simulation.py`
- interpolacao de politicas em `zhang2005/interpolation.py`

Caminhos naturais:

1. Numba nas rotinas de VFI e simulacao.
2. F2PY envolvendo os Fortran originais.
3. Rodadas paralelas para simulacoes independentes.

## Problemas Comuns

### `python: can't open file 'Untitled-1'`

Isso significa que o VS Code tentou rodar uma aba nao salva. Salve o arquivo ou rode diretamente um script existente, por exemplo:

```bash
.venv/bin/python scripts/run_calibration.py
```

### `ModuleNotFoundError: No module named numpy`

O ambiente virtual nao esta instalado ou nao esta sendo usado. Rode:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

E use sempre:

```bash
.venv/bin/python ...
```

### Figuras nao aparecem

Os scripts salvam arquivos em disco, nao abrem janela grafica. Veja:

```text
outputs/figures/png/
outputs/figures/pdf/
```

### Tabelas com momentos extremos

A simulacao reduzida pode gerar algumas observacoes de valor de firma muito proximas de zero. As tabelas model-only filtram precos abaixo de `1e-6` nas ordenacoes de book-to-market para evitar que uma firma quase-zero domine os deciles. Isso e uma decisao pratica para a fase reduzida; no benchmark final a qualidade da simulacao deve ser reavaliada.

## Proximos Passos

1. Criar uma rotina de benchmark medio:
   - `N = 500`
   - `periods = 421`
   - `ks_iterations = 10`
   - `vfi_iterations >= 250`

2. Implementar aceleracao com Numba.

3. Rodar Table IV de estatica comparativa:
   - baixa/alta volatilidade idiossincratica;
   - ajuste rapido/lento;
   - baixo/alto custo fixo;
   - simetria nos custos de ajuste;
   - preco de risco constante.

4. Refinar Figura 4 estrutural:
   - calcular retorno esperado por estado;
   - comparar baixo `z` versus alto `z`;
   - plotar contra produtividade agregada `x`.

5. Criar notebooks narrativos finais.

## Referencias Locais

- Artigo em PDF:
  - `The Journal of Finance - 2005 - ZHANG - The Value Premium.pdf`

- Codigo original:
  - `vpCode/`

- Notas tecnicas do port:
  - `docs/zhang2005_porting_notes.md`
