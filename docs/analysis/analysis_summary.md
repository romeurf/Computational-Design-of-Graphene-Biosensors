# Análise exploratória — pipeline GFET
Fonte: `FINAL_PROBES_ALL.csv` (8455 probes; 8455 próprias, 0 Beatriz/literatura).

## 1. seqfold — distribuição de ΔG MFE e limiar
- Probes consideradas (próprias, pass_basic, ΔG definido): **3064**  (15 com ΔG implausível >+50 kcal/mol = artefacto seqfold, excluídos das figuras/estatísticas)
- ΔG MFE (plausível): min **-5.40**, mediana **-0.20** kcal/mol
- Percentis de ΔG (cauda mais estável = mais negativa): P1=-3.20, P5=-2.60, P10=-2.10, P25=-1.10, P50=-0.20

- Aprovação a vários limiares (ΔG ≥ limiar):

| limiar | passam (global) | % | nuc | rmpM | lytA | oprL | algD | frdB |
|---|---|---|---|---|---|---|---|---|
| 0 | 1470 | 48.0 | 82 | 413 | 417 | 86 | 155 | 317 |
| -1 | 2265 | 73.9 | 114 | 630 | 563 | 233 | 269 | 456 |
| -2 | 2745 | 89.6 | 119 | 762 | 726 | 289 | 348 | 501 |
| -3 | 3025 | 98.7 | 143 | 841 | 803 | 322 | 400 | 516 |
| -4 | 3055 | 99.7 | 143 | 841 | 807 | 342 | 403 | 519 |
| -5 | 3059 | 99.8 | 143 | 841 | 807 | 344 | 403 | 521 |
| -6 | 3064 | 100.0 | 143 | 841 | 807 | 349 | 403 | 521 |

**Recomendação (para discussão):** limiar atual -2.6 kcal/mol → 96.1% passam (população pass_basic, que já tende a ter pouca estrutura). O P5 da distribuição ≈ **-2.6 kcal/mol** marca as ~5% mais estruturadas; o limiar atual está alinhado com o P5 — defensável estatisticamente. Ajustar conforme discussão com o orientador.

## 2. Tamanhos & hits por gene + normalidade

| gene | nº seqs (hits) | comp. min–máx | média±DP | Shapiro p | normal? |
|---|---|---|---|---|---|
| nuc | 77 | 203–333 | 245.5±23.3 | 0.0 | não |
| rmpM | 76 | 870–1185 | 1089.2±110.7 | 0.0 | não |
| lytA | 87 | 700–1132 | 902.2±113.8 | 0.0003 | não |
| oprL | 24 | 413–562 | 481.8±31.2 | 0.043 | não |
| algD | 12 | 526–657 | 563.7±32.0 | 0.0007 | não |
| frdB | 82 | 489–489 | 489.0±0.0 | n/a | n/a |

_Shapiro-Wilk: p < 0.05 ⇒ rejeita normalidade. n/a quando n<3 ou variância nula._

- Comprimento das probes próprias: n=8455, intervalo 18–28 nt, média 22.9 (discreto/limitado 18–28 nt → não-normal por construção).

## 3. Diversidade entre sequências recuperadas (sem alinhamento, k-mer)
- Distância = 1 − similaridade do cosseno entre vetores de frequência de 4-mers.

| gene | nº seqs | diversidade média (0=idênticas, →1 diversas) |
|---|---|---|
| nuc | 77 | 0.0594 |
| rmpM | 76 | 0.0793 |
| lytA | 87 | 0.1868 |
| oprL | 24 | 0.0503 |
| algD | 12 | 0.0966 |
| frdB | 82 | 0.0254 |

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
- 8455 probes descritas → `output/analysis/probe_descriptors.csv` (29 colunas: composição, GC/AT, purina/pirimidina, entropia, 16 dinucleótidos).
- _Descritores de estruturas 3D: adiados (sem estruturas em disco) — ficam como próximo passo._

## 6. Parâmetros por gene (revisão / transparência)
- **Auto por gene:** o comprimento é selecionado automaticamente (cluster dominante ±25% pós-fetch) → genes novos não precisam de afinar min_len/max_len à mão.
- **Fixos (decisões biológicas, override por gene em TARGETS):** cons_min, GC, Tm.

| gene | n seqs | comp. usado | filtro min–max | tol | cons_min | GC | Tm_min |
|---|---|---|---|---|---|---|---|
| nuc | 77 | 203–333 | 200–3000 | ±0.25 | 0.85 | 0.38–0.60 | 52.0 |
| rmpM | 76 | 870–1185 | 100–1200 | ±0.25 | 0.8 | 0.40–0.65 | 53.0 |
| lytA | 87 | 700–1132 | 700–1300 | ±0.25 | 0.7 | 0.40–0.60 | 53.0 |
| oprL | 24 | 413–562 | 300–2000 | ±0.25 | 0.85 | 0.40–0.70 | 53.0 |
| algD | 12 | 526–657 | 500–2500 | ±0.25 | 0.85 | 0.40–0.70 | 53.0 |
| frdB | 82 | 489–489 | 200–2000 | ±0.25 | 0.85 | 0.38–0.60 | 53.0 |

_O comprimento usado é o cluster dominante já filtrado pelo pipeline. cons_min/GC/Tm permanecem defaults informados (ex.: lytA cons 0.70 por diversidade alélica — Whatmore 2000; oprL/algD GC≤0.70 — Stover 2000)._