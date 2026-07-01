# Resumo — Duas adições e resultados

**Romeu Fernandes (PG45861) · Universidade do Minho · 2026-06-16**

---

## 1. Integração de probes externas (IPLEX) no pipeline

**O que fiz:** o pipeline passou a aceitar um ficheiro externo (xlsx/csv com as colunas
`job_name, sequence, species`) e a **pontuar essas probes com as mesmas métricas** das
minhas (Tm, GC, hairpin, homodímero, seqfold, No-fold), entrando na **mesma tabela
consolidada** — sem caminho separado (`--with-reference`). A conservação (PPI) fica
explicitamente **não aplicável**, porque só é definível para probes desenhadas a partir de
um alinhamento; probes já feitas não têm região de origem.

**Resultado — as 74 probes de `data/sequencias_iplex.xlsx` no meu pipeline:**

Funil de triagem (critérios globais):

| Passo | Passam |
|---|---|
| Total | 74 |
| Tm 53–72 °C | 56 |
| + GC 40–60% | 36 |
| + hairpin > −2 kcal/mol | 35 |
| + homodímero > −5 kcal/mol | **29** (passam a triagem básica = 39%) |
| + seqfold ≥ −2,6 kcal/mol | **29** (passam o crivo completo) |

**29 das 74 (39%)** passam o crivo completo. A maior perda é no **GC** (20 falham a janela
40–60%). **Nenhuma** é rejeitada pelo seqfold — as que passam já têm pouca estrutura secundária.

Comparação (médias) com as minhas probes aprovadas:

| | IPLEX | Minhas |
|---|---|---|
| Tm (°C) | 56,3 | 59,3 |
| GC | 0,50 | 0,52 |
| No-fold | 80 | 83 |

→ as IPLEX ficam ligeiramente abaixo (Tm e No-fold mais baixos), coerente com serem um painel
mais antigo e heterogéneo (inclui alvos virais/fúngicos) que não foi otimizado pelo meu crivo.

**Como aparecem no pipeline:** ao correr `python pipeline.py --max-seqs 100 --with-reference`,
estas 74 probes são pontuadas e entram no `output/FINAL_PROBES_ALL.csv` com **`fonte = referencia`**,
lado a lado com as minhas e com as mesmas colunas → permite **comparação direta na mesma tabela**.

## 2. Novas métricas de diversidade (sem alinhamento)

Todas são **baseadas em k-mers, sem alinhamento** (o objetivo era avaliar a diversidade
*sem alinhar* as sequências), pelo que funcionam mesmo com comprimentos variáveis e
sequências divergentes.

### 2.1 Porquê estas métricas (e não outras)
São **complementares** — cada uma mede um aspeto diferente e, em conjunto, dão uma leitura
robusta que uma média sozinha não daria:
- **Cosseno** pesa as *frequências* dos k-mers; **Jaccard** usa só *presença/ausência*
  (o "vocabulário" de k-mers). Se as duas concordam, a diversidade é real, não um artefacto
  de uma métrica.
- **Desvio-padrão** e **máximo** do cosseno mostram a *forma* da distribuição (subgrupos/outliers
  que a média esconde).
- **% únicas** deteta o artefacto dos *duplicados* (o NCBI devolve muitas cópias).

Preteri a distância **Euclidiana** (sensível ao comprimento das sequências) e as medidas
**baseadas em alinhamento** (π, identidade média — violam o requisito "sem alinhamento",
são lentas e falham em conjuntos divergentes).

### 2.2 Como ler cada métrica

| Métrica | O que mede | Valor **alto** → | Valor **baixo** → |
|---|---|---|---|
| Distância de cosseno média | distância composicional média (frequências de k-mers) | conjunto **diverso** | conjunto **conservado/semelhante** |
| Desvio-padrão da distância de cosseno | dispersão da diversidade entre pares | **heterogéneo** (há subgrupos/outliers) | diversidade **uniforme** |
| Distância de cosseno máxima | o par de sequências mais divergente | existe pelo menos um **outlier**/subgrupo distinto | até o par mais diferente é **próximo** |
| Distância de Jaccard média | distância por presença/ausência de k-mers | partilham **poucos** k-mers (vocabulário diferente) | partilham **quase todos** os k-mers |
| Sequências únicas (%) | redundância (duplicados exatos) | sequências **todas distintas** | **muitos duplicados** no conjunto |

### 2.3 Resultados (por gene)

| gene | nº seqs | cosseno médio | desvio-padrão | cosseno máximo | Jaccard médio | sequências únicas (%) |
|---|---|---|---|---|---|---|
| nuc  | 77 | 0,059 | 0,074 | 0,425 | 0,163 | 61,0 |
| rmpM | 76 | 0,079 | 0,076 | 0,247 | 0,069 | 50,0 |
| lytA | 87 | **0,187** | 0,144 | 0,440 | 0,125 | 93,1 |
| oprL | 24 | 0,050 | 0,084 | 0,358 | 0,096 | 91,7 |
| algD | 12 | 0,097 | 0,140 | 0,332 | 0,117 | **100,0** |
| frdB | 82 | 0,025 | 0,013 | 0,052 | 0,061 | **35,4** |

**Conclusões:**
- **lytA** é o mais diverso (cosseno médio 0,187) — coerente com a sua diversidade alélica.
- **frdB** é o mais conservado e redundante (cosseno 0,025 e só **35% de sequências únicas**
  → muitos duplicados no NCBI).
- **algD** tem 100% de sequências únicas (poucas, mas todas distintas) → explica a sua
  instabilidade e reforça a necessidade da seleção automática de comprimento.
- **Cosseno e Jaccard concordam** entre si → leitura robusta da diversidade, sem alinhamento.
