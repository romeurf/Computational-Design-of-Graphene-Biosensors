# Justificação do Pipeline de Design de Probes para Biossensores GFET

**Autor:** Romeu Fernandes (PG45861) · Universidade do Minho
**Data:** 2026-06-16
**Código:** `pipeline.py` · análise: `scripts/analysis.py` · 3D: `colab_boltz2_batch.ipynb`

---

## 1. Objetivo

Desenhar *in silico* probes de ssDNA para deteção de 6 patógenos bacterianos em
biossensores GFET (Graphene Field-Effect Transistor). O pipeline parte de sequências
genómicas públicas e filtra, de forma fundamentada e reprodutível, até um conjunto
reduzido de candidatas ótimas, validadas por previsão de estrutura 3D, minimizando o
tempo de cálculo gasto em más candidatas.

## 2. Alvos

| Gene | Organismo | Grupo | Função |
|---|---|---|---|
| nuc | *Staphylococcus aureus* | A | termonuclease |
| rmpM | *Neisseria meningitidis* | A | proteína de membrana externa (classe 4) |
| lytA | *Streptococcus pneumoniae* | B | autolisina |
| oprL | *Pseudomonas aeruginosa* | B | lipoproteína associada ao peptidoglicano |
| algD | *Pseudomonas aeruginosa* | B | GDP-manose 6-desidrogenase |
| frdB | *Haemophilus influenzae* | B | subunidade Fe-S da fumarato redutase |

## 3. Etapas do pipeline e justificação

1. **Recuperação NCBI** — `esearch`/`efetch` por gene+organismo, nº configurável
   (`--max-seqs`, default usado: 100). *Porquê 100:* ver §5 (rarefação).
2. **Cluster de comprimento automático** — após o fetch, mantém-se só o cluster de
   comprimento dominante (±25% em torno do comprimento mais denso). *Porquê:* misturar
   fragmentos parciais com registos completos/sobredimensionados gera alinhamentos com
   muitos gaps e poucas janelas conservadas. A seleção é **data-driven** (não há
   `min_len`/`max_len` afinados à mão por gene), pelo que funciona para genes novos sem
   ajuste. *Impacto medido:* resolveu o algD (0 → ~400 probes) e melhorou o nuc.
3. **Alinhamento múltiplo (MAFFT)** — `--auto` (FFT-NS-2); modo sequência única quando
   só há 1 sequência.
4. **Conservação = PPI** (Percentage of Pairwise Identity) — portado do **ViruScope**
   (Ana Lima, `anasfplima/ViruScope`). É a identidade par-a-par média por coluna do
   alinhamento, com pesos IUPAC para bases ambíguas — métrica mais rigorosa que a simples
   frequência da base mais comum. **PPI3 (identidade dos 3 nt terminais) foi excluída de
   propósito**: ao contrário de primers, as probes não são estendidas por polimerase,
   pelo que a identidade da extremidade 3′ é irrelevante.
5. **Triagem termoquímica (primer3-py)** — Tm (nearest-neighbour), GC%, hairpin ΔG,
   homodímero ΔG. Limiares na §4.
6. **Estrutura secundária (seqfold)** — ΔG MFE a 37 °C; rejeita probes com estrutura
   demasiado estável (que hibridizariam mal). Limiar recalibrado (§4).
7. **No-fold score** (ViruScope) — transformação logística do ΔG (hairpin/homodímero/
   seqfold) num score contínuo 0–100, alternativa ao pass/fail rígido.
8. **Seleção para 3D** — top-N por gene **por qualidade** = média de PPI e No-fold
   (não por Tm). *Porquê:* os passes já filtram Tm/GC/estrutura; ordenar os sobreviventes
   pelos sinais de qualidade do pipeline gasta o tempo de GPU nas melhores candidatas.
9. **Previsão 3D (Boltz-2, Colab GPU)** — confidence, pTM, pLDDT por probe + **matriz PAE**
   (estilo AlphaFold3). *Nota:* pTM/PAE são calibrados para proteínas; em ssDNA curto o
   confidence e o pLDDT são os indicadores úteis.
10. **Shortlist final** — fusão dos resultados Boltz com os metadados → ranking combinando
    qualidade de sequência (PPI+No-fold) e 3D (confidence/pLDDT).

## 4. Parâmetros fixos e fundamentação na literatura

| Parâmetro | Valor | Fundamento |
|---|---|---|
| Comprimento da probe | 18–28 nt | Wetmur 1991; GFET usam tipicamente 20-mer ssDNA |
| Tm | 53–72 °C (52 para AT-rich) | SantaLucia & Hicks 2004; deriva auto da T do ensaio via `--assay-temp` (T+15…+35) |
| GC | 40–60% (38–60 AT-rich; ≤70 P. aeruginosa) | janela ótima 40–60% (PremierBiosoft/UNLV); P. aeruginosa ~67% GC genómico (Stover 2000) |
| Hairpin ΔG | > −2.0 kcal/mol | guias IDT/primer3 (−1.5 a −3) |
| Homodímero ΔG | > −5.0 kcal/mol | guias IDT (self-dimer ~−6) |
| Conservação (PPI) | ≥ 0.85 (rmpM 0.80; lytA 0.70) | regiões conservadas para deteção pan-estirpe; lytA relaxado por diversidade alélica (Whatmore 2000) |
| seqfold ΔG MFE | ≥ −2.6 kcal/mol | **P5 da distribuição observada** (fixo/reprodutível); coincide com a guia "hairpin −3"; recalibrado de −6.0 (que deixava passar 100%) |

*Decisão:* Tm/GC/conservação são **defaults biológicos fixos** (com override por gene
quando justificado) — **não** são auto-derivados dos dados, para não baixar a fasquia em
genes mal-conservados. Apenas o **comprimento** é data-driven (§3.2), por ser uma questão
de qualidade do alinhamento e não um requisito biológico.

## 5. Como o N (nº de sequências) foi decidido — rigor, não arbitrário

Análise de **rarefação** (`scripts/analysis.py`): subamostrando N sequências por gene e
medindo a diversidade k-mer e a riqueza (k-mers distintos), estas **saturam por volta de
N≈50–75** no gene mais exigente (nuc); a maioria satura a 5–20. Logo **N≈100 por gene é
suficiente** (com margem) — mais do que isso quase não acrescenta informação nova. O fator limitante
real são os genes com poucas sequências no NCBI (algD ~19, oprL ~24), não o limite N.
*Exceção:* `algD` tem override `max_seqs=200` (usar todas as disponíveis; instável a 100).

## 6. Decisões de design registadas

- **Conservação = só PPI** (sem PPI3) — probes não são estendidas por polimerase.
- **Comprimento automático** por cluster dominante — sem afinação manual por gene.
- **Limiares recalibrados a partir dos dados/literatura** (gc_max 0.60; seqfold −2.6),
  mantidos **fixos** (reprodutíveis) em vez de recalculados a cada corrida.
- **Janela Tm derivável** da temperatura do ensaio (`--assay-temp`); por defeito 53–72 °C.
- **Seleção para Boltz por qualidade** (PPI+No-fold), não por Tm.
- **Merge de probes de referência é opcional** (`--with-reference`) — o core corre sobre as probes próprias.

## 7. Resultados (corrida N=100, 2026-06-16)

Funil global: **8455 janelas candidatas → 3078 passam triagem básica → 2943 passam
seqfold → 30 selecionadas para Boltz** (top 5/gene por qualidade).

Validação 3D (Boltz-2): **4 GOOD, 11 MODERATE, 15 LOW** (confidence). Melhores candidatas
(fortes em sequência *e* em 3D):

| probe | gene | PPI | No-fold | confidence | pLDDT |
|---|---|---|---|---|---|
| p5972 | algD | 0.93 | 98 | 0.731 | 0.901 |
| p7716 | frdB | 0.99 | 100 | 0.723 | 0.822 |
| p5938 | algD | 0.90 | 96 | 0.721 | 0.887 |
| p8448 | frdB | 1.00 | 98 | 0.716 | 0.845 |

Shortlist completa: `output/colab_boltz2/boltz2_shortlist_ranked.csv`.
Análise exploratória (distribuições, normalidade, diversidade): `output/analysis/` + `output/figures/`.

## 8. Reprodutibilidade — como correr

```bash
python pipeline.py --max-seqs 100 --colab 5      # pipeline completo → boltz2_inputs.zip
# (Colab) correr colab_boltz2_batch.ipynb com o boltz2_inputs.zip → resultados
python pipeline.py --merge-boltz output/colab_boltz2/boltz2_results_summary.csv  # shortlist
python scripts/analysis.py                       # figuras + análise exploratória
```

*Nota de reprodutibilidade:* o NCBI não devolve sempre o mesmo conjunto, pelo que as
contagens por gene podem variar ±poucos % entre corridas. Para reprodutibilidade estrita,
fixar os accession IDs.

## 9. Referências

- SantaLucia & Hicks 2004 — termodinâmica nearest-neighbour (Tm).
- Wetmur 1991 — comprimento de probes (18–28 nt).
- Zadeh et al. 2011 (NUPACK) — energia de estrutura secundária.
- IDT OligoAnalyzer — GC, hairpin/dímero ΔG.
- Stover et al. 2000 — conteúdo GC genómico de *P. aeruginosa*.
- Whatmore et al. 2000 — diversidade alélica de *lytA* em *S. pneumoniae*.
- Lima et al. — ViruScope (`github.com/anasfplima/ViruScope`): PPI, No-fold score.
- Wohlwend et al. 2024 — Boltz-2 (previsão de estrutura 3D).
