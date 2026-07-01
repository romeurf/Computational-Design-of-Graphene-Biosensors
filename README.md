# Computational Design of Graphene Biosensors

![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python&logoColor=white)
![Status](https://img.shields.io/badge/Status-In%20Development-yellow)
![University](https://img.shields.io/badge/University-Minho-darkgreen)

> Master's Thesis — Bioinformatics / Computational Biology  
> University of Minho · PG45861  
> Computational pipeline for the *in silico* design of DNA probes for **Graphene Field-Effect Transistor (GFET) biosensors** targeting bacterial pathogens.

---

## Bacterial Targets

| Gene | Organism | Group |
|------|----------|-------|
| `nuc` | *Staphylococcus aureus* | A |
| `rmpM` | *Neisseria meningitidis* | A |
| `lytA` | *Streptococcus pneumoniae* | B |
| `oprL` | *Pseudomonas aeruginosa* | B |
| `algD` | *Pseudomonas aeruginosa* | B |
| `frdB` | *Haemophilus influenzae* | B |

---

## Pipeline

```
NCBI Entrez
    │
    ▼
MAFFT (múltiplo alinhamento)
    │
    ▼
Janelas candidatas (conservação ≥ 85%, gap ≤ 20%, 18–28 nt)
    │
    ▼
Scoring básico — primer3-py (Tm, GC, hairpin ΔG, homodimer ΔG)
    │
    ▼
seqfold — estrutura secundária MFE (substituto open-source do NUPACK)
    │
    ▼
Boltz-2 — predição estrutura 3D + RMSD entre réplicas
    │
    ▼
output/FINAL_PROBES_ALL.csv   ← probes próprias (+ referência opcional) consolidadas
```

---

## Critérios por gene

| Parâmetro | Global | nuc/frdB (AT-rico) | oprL/algD (GC-rico) |
|-----------|--------|--------------------|---------------------|
| Tm mín | 53 °C | 52 °C | 53 °C |
| GC | 40–65% | 38–60% | 40–70% |
| Hairpin ΔG | > −2 kcal/mol | = | = |
| Homodimer ΔG | > −5 kcal/mol | = | = |
| seqfold MFE | > −3 kcal/mol | = | = |
| Bases emparelhadas | < 15% | = | = |

Ref: SantaLucia & Hicks 2004 · Wetmur 1991 · IDT OligoAnalyzer · Stover 2000

---

## Estrutura do repositório

```
/
├── pipeline.py              ← entrada principal
├── scripts/
│   ├── boltz2_predict.py    ← predição 3D Boltz-2 (standalone)
│   └── nucleofold_predict.py← wrapper NucleoFold3D (requer WSL/Linux)
├── data/
│   ├── reference_probes_IPLEXMED.xlsx
│   └── tabela_targets_probes_updated.xlsx
├── docs/
│   ├── Protocolo_Fase0_Alinhamento.docx
│   ├── The_GFET_Computational_Blueprint.pdf
│   ├── Thesis_Plan_PG45861.pdf
│   └── documentacao_tecnica_fetch_and_align.docx
└── output/                  ← gerado em runtime (gitignored)
    ├── alignments/
    ├── structures/
    └── FINAL_PROBES_ALL.csv
```

---

## Instalação

```bash
pip install biopython primer3-py pandas requests pyyaml seqfold openpyxl
```

MAFFT em PATH (instalar de https://mafft.cbrc.jp/alignment/software/).

Boltz-2 (opcional, para predição 3D):
```bash
pip install boltz
```

---

## Correr o pipeline

```bash
# Pipeline completo (todos os targets)
python pipeline.py

# Sem predição 3D (mais rápido)
python pipeline.py --no-3d

# Sem seqfold (só critérios básicos)
python pipeline.py --no-nupack --no-3d

# Target específico
# editar TARGETS em pipeline.py para comentar os que não interessam
```

### Outputs gerados

| Ficheiro | Descrição |
|---------|-----------|
| `output/FINAL_PROBES_ALL.csv` | Todas as probes (próprias + referência) com scores |
| `output/alignments/<gene>/<gene>_probes_scored.tsv` | Probes por gene com todas as métricas |
| `output/alignments/<gene>/<gene>_viroscope_probes.fasta` | FASTA das probes PASS |
| `output/alignments/<gene>/aligned.fasta` | Alinhamento múltiplo (Geneious/AliView) |
| `output/structures/` | Estruturas 3D (PDB/CIF) |

### Predição 3D standalone (após pipeline)

```bash
# Boltz-2 (Windows nativo)
python scripts/boltz2_predict.py --filter nupack --samples 3

# NucleoFold3D (requer WSL/Linux com nsp + AmberTools)
python scripts/nucleofold_predict.py --prepare-only   # prepara CSV
# depois em WSL: python NucleoFold3D.py --csv output/nucleofold_input.csv
python scripts/nucleofold_predict.py --collect        # recolhe resultados
```

---

## Referências dos limiares

| Critério | Fonte |
|---------|-------|
| Tm nearest-neighbour | SantaLucia & Hicks 2004 |
| Comprimento probe 18–28 nt | Wetmur 1991 |
| seqfold MFE (proxy NUPACK) | Zadeh et al. 2011 (NUPACK) |
| GC 40–65%, homodimer > −5 | IDT OligoAnalyzer |
| GC P. aeruginosa até 70% | Stover et al. 2000 (PAO1) |
| Boltz-2 estrutura 3D | Wohlwend et al. 2024 |
| RMSD Biopython | Hamelryck & Manderick 2003 |

---

*Academic use only — University of Minho 2025–2026*
