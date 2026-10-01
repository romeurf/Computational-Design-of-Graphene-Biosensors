# Análise exploratória — pipeline GFET
Fonte: `FINAL_PROBES_ALL.csv` (7735 probes; 7661 próprias, 74 referência/literatura).

## 1. seqfold — distribuição de ΔG MFE e limiar
- Probes consideradas (próprias, pass_basic, ΔG definido): **2554**  (26 com ΔG implausível >+50 kcal/mol = artefacto seqfold, excluídos das figuras/estatísticas)
- ΔG MFE (plausível): min **-5.10**, mediana **0.10** kcal/mol
- Percentis de ΔG (cauda mais estável = mais negativa): P1=-2.90, P5=-2.10, P10=-1.50, P25=-0.70, P50=0.10

- Aprovação a vários limiares (ΔG ≥ limiar):

| limiar | passam (global) | % | nuc | rmpM | lytA | oprL | algD | frdB |
|---|---|---|---|---|---|---|---|---|
| 0 | 1485 | 58.1 | 119 | 172 | 539 | 89 | 131 | 435 |
| -1 | 2148 | 84.1 | 158 | 295 | 669 | 225 | 217 | 584 |
| -2 | 2406 | 94.2 | 164 | 320 | 778 | 254 | 276 | 614 |
| -3 | 2541 | 99.5 | 181 | 322 | 851 | 254 | 308 | 625 |
| -4 | 2551 | 99.9 | 181 | 322 | 856 | 258 | 308 | 626 |
| -5 | 2551 | 99.9 | 181 | 322 | 856 | 258 | 308 | 626 |
| -6 | 2554 | 100.0 | 181 | 322 | 856 | 261 | 308 | 626 |

**Recomendação (para discussão):** limiar atual -2.1 kcal/mol → 95.5% passam (população pass_basic, que já tende a ter pouca estrutura). O P5 da distribuição ≈ **-2.1 kcal/mol** marca as ~5% mais estruturadas; o limiar atual está alinhado com o P5 — defensável estatisticamente. Ajustar conforme discussão com o orientador.

## 2. Tamanhos & hits por gene + normalidade

| gene | nº seqs (hits) | comp. min–máx (bp) | range (bp) | média±DP | Shapiro p | normal? |
|---|---|---|---|---|---|---|
| nuc | 77 | 203–333 | 130 | 245.5±23.3 | 0.0 | não |
| rmpM | 76 | 870–1185 | 315 | 1089.2±110.7 | 0.0 | não |
| lytA | 87 | 700–1132 | 432 | 902.2±113.8 | 0.0003 | não |
| oprL | 24 | 413–562 | 149 | 481.8±31.2 | 0.043 | não |
| algD | 12 | 526–657 | 131 | 563.7±32.0 | 0.0007 | não |
| frdB | 82 | 489–489 | 0 | 489.0±0.0 | n/a | n/a |

_Shapiro-Wilk: p < 0.05 ⇒ rejeita normalidade. n/a quando n<3 ou variância nula._

- Comprimento das probes próprias: n=7661, intervalo 18–28 nt, média 22.9 (discreto/limitado 18–28 nt → não-normal por construção).

## 3. Diversidade entre sequências recuperadas (sem alinhamento, k-mer)
- Vetores de frequência de 4-mers por sequência; medidas complementares (todas **sem alinhamento**): **cosseno** (1−similaridade — média/DP/máx), **Jaccard** (presença/ausência de k-mers) e **% de sequências únicas**.

| gene | nº seqs | cosseno médio | cosseno DP | cosseno máx | Jaccard médio | % únicas |
|---|---|---|---|---|---|---|
| nuc | 77 | 0.0594 | 0.0735 | 0.4246 | 0.1626 | 61.0 |
| rmpM | 76 | 0.0793 | 0.0763 | 0.2469 | 0.0693 | 50.0 |
| lytA | 87 | 0.1868 | 0.1443 | 0.4395 | 0.1254 | 93.1 |
| oprL | 24 | 0.0503 | 0.0841 | 0.3578 | 0.096 | 91.7 |
| algD | 12 | 0.0966 | 0.1398 | 0.3318 | 0.1169 | 100.0 |
| frdB | 82 | 0.0254 | 0.0132 | 0.0517 | 0.061 | 35.4 |
_0 = idênticas, →1 = diversas. DP/máx mostram a dispersão da diversidade; % únicas indica redundância (duplicados) no conjunto recuperado._

## 4. Rarefação — quantas sequências usar (escalar com rigor)
- Para cada gene, subamostram-se N sequências e mede-se a diversidade k-mer média e a riqueza (k-mers distintos). Quando saturam, mais sequências acrescentam pouco.

| gene | nº disponível | N de saturação (riqueza ≥95% do máx) |
|---|---|---|
| nuc | 77 | 50 |
| rmpM | 76 | 5 |
| lytA | 87 | 5 |
| oprL | 24 | 20 |
| algD | 12 | 10 |
| frdB | 82 | 5 |

**Recomendação de N:** a diversidade/riqueza satura por volta de **N ≈ 50** sequências para o gene mais exigente. Usar N ≈ 50–100 por gene é suficiente (mais do que isso acrescenta pouca informação nova). Genes com poucas sequências no NCBI (ex.: algD) são o fator limitante real, não o cap.

## 5. Descritores de sequência (alternativa ao PyBioMed)
- 7735 probes descritas → `output/analysis/probe_descriptors.csv` (29 colunas: composição, GC/AT, purina/pirimidina, entropia, 16 dinucleótidos).
- _Descritores de estruturas 3D: adiados (sem estruturas em disco) — ficam como próximo passo._

## 6. Parâmetros por gene (revisão / transparência)
- **Auto por gene:** o comprimento é selecionado automaticamente (cluster dominante ±25% pós-fetch) → genes novos não precisam de afinar min_len/max_len à mão.
- **Fixos (decisões biológicas, override por gene em TARGETS):** cons_min, GC, Tm.

| gene | n seqs | comp. usado | filtro min–max | tol | cons_min | GC | Tm_min |
|---|---|---|---|---|---|---|---|
| nuc | 77 | 203–333 | 200–3000 | ±0.25 | 0.85 | 0.38–0.60 | 52.0 |
| rmpM | 76 | 870–1185 | 100–1200 | ±0.25 | 0.85 | 0.40–0.60 | 53.0 |
| lytA | 87 | 700–1132 | 700–1300 | ±0.25 | 0.7 | 0.40–0.60 | 53.0 |
| oprL | 24 | 413–562 | 300–2000 | ±0.25 | 0.85 | 0.40–0.70 | 53.0 |
| algD | 12 | 526–657 | 500–2500 | ±0.25 | 0.85 | 0.40–0.70 | 53.0 |
| frdB | 82 | 489–489 | 200–2000 | ±0.25 | 0.85 | 0.38–0.60 | 53.0 |

_O comprimento usado é o cluster dominante já filtrado pelo pipeline. cons_min/GC/Tm permanecem defaults informados (ex.: lytA cons 0.70 por diversidade alélica — Whatmore 2000; oprL/algD GC≤0.70 — Stover 2000)._