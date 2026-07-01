# Resumo — Duas adições e resultados

**Romeu Fernandes (PG45861) · Universidade do Minho · 2026-06-16**

---

## 1. Integração de probes externas (IPLEX) no pipeline

**O que fiz:** o pipeline passou a aceitar um ficheiro externo (xlsx/csv com as colunas
`job_name, sequence, species`) e a **pontuar essas probes com as mesmas métricas** das
minhas (Tm, GC, hairpin, homodímero, seqfold, No-fold), entrando na **mesma tabela
consolidada** — sem caminho separado (`--with-reference`). A conservação (PPI) fica
explicitamente **N/A**, porque só é definível para probes desenhadas a partir de um
alinhamento; probes já feitas não têm região de origem.

**Resultado (74 probes IPLEX):**
- **29 / 74 (39%)** passam a triagem básica com os critérios globais.
- Comparação (médias) com as minhas probes aprovadas:

  | | IPLEX | Minhas |
  |---|---|---|
  | Tm (°C) | 56,3 | 59,3 |
  | GC | 0,50 | 0,52 |
  | No-fold | 80,2 | 82,9 |

  → as IPLEX ficam ligeiramente abaixo (Tm e No-fold mais baixos), coerente com serem um
  painel mais antigo/heterogéneo (inclui alvos virais/fúngicos).

## 2. Novas métricas de diversidade (sem alinhamento)

**O que fiz:** além da diversidade média por cosseno, acrescentei **desvio-padrão e máximo
do cosseno**, **distância de Jaccard** (presença/ausência de k-mers) e **% de sequências
únicas** — todas sem alinhamento.

**Resultado (por gene):**

| gene | nº seqs | cosseno médio | cosseno DP | cosseno máx | Jaccard médio | % únicas |
|---|---|---|---|---|---|---|
| nuc  | 77 | 0,059 | 0,074 | 0,425 | 0,163 | 61,0 |
| rmpM | 76 | 0,079 | 0,076 | 0,247 | 0,069 | 50,0 |
| lytA | 87 | **0,187** | 0,144 | 0,440 | 0,125 | 93,1 |
| oprL | 24 | 0,050 | 0,084 | 0,358 | 0,096 | 91,7 |
| algD | 12 | 0,097 | 0,140 | 0,332 | 0,117 | **100,0** |
| frdB | 82 | 0,025 | 0,013 | 0,052 | 0,061 | **35,4** |

**Conclusões:**
- **lytA é o gene mais diverso** (cosseno médio 0,187) — coerente com a sua diversidade alélica.
- **frdB é o mais conservado e redundante** (cosseno 0,025 e só **35% de sequências únicas**
  → muitos duplicados no NCBI).
- **algD** tem 100% de sequências únicas (poucas mas todas distintas), o que explica a sua
  instabilidade — reforça a necessidade da seleção automática de comprimento.
- As várias medidas concordam entre si (cosseno ↔ Jaccard), dando uma leitura robusta da
  diversidade sem precisar de alinhamento.
