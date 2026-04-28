# 🧬 Computational Design of Graphene Biosensors

![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)
![License](https://img.shields.io/badge/License-Academic%20Use%20Only-lightgrey)
![Status](https://img.shields.io/badge/Status-In%20Development-yellow)
![University](https://img.shields.io/badge/University-Minho-darkgreen?logo=academia&logoColor=white)
![Thesis](https://img.shields.io/badge/Thesis-Bioinformatics%20%2F%20CompBio-blueviolet)
![MAFFT](https://img.shields.io/badge/Alignment-MAFFT%207.526-orange)
![primer3](https://img.shields.io/badge/Thermodynamics-primer3--py-red)
![NCBI](https://img.shields.io/badge/Data-NCBI%20Entrez-0070c0?logo=databricks&logoColor=white)

> **Master's Thesis — Bioinformatics / Computational Biology**  
> University of Minho · PG45861  
> Computational pipeline for the *in silico* design of DNA probes for **Graphene Field-Effect Transistor (GFET) biosensors** targeting bacterial pathogens.

---

## 🔬 Overview

This project implements a full computational workflow to design and score oligonucleotide probes for GFET biosensors capable of detecting clinically relevant bacterial pathogens. The pipeline integrates sequence acquisition, multiple sequence alignment, conservation analysis, and thermodynamic scoring — all in Python.

The work is structured in two main phases:

| Phase | Description |
|-------|-------------|
| **Phase 0** | Sequence acquisition (NCBI) → MAFFT alignment → conservation analysis → probe candidate extraction → ViruScope-style thermodynamic scoring |
| **Phase 1** | 3D DNA modelling (3dDNA) → molecular dynamics (AMBER) → molecular docking |

---

## 🦠 Bacterial Targets

Six clinically validated gene markers across two groups:

| Gene | Organism | Group | Clinical Relevance |
|------|----------|-------|-------------------|
| `nuc` | *Staphylococcus aureus* | A | Thermostable nuclease — GFET benchmark (Purwidyantri 2021) |
| `rmpM` | *Neisseria meningitidis* | A | Class 4 OMP — electrochemical benchmark (Appaturi 2013) |
| `lytA` | *Streptococcus pneumoniae* | B | Autolysin LytA — gold-standard clinical marker |
| `oprL` | *Pseudomonas aeruginosa* | B | Lipoprotein OprL — gold-standard |
| `algD` | *Pseudomonas aeruginosa* | B | Mucoidal/biofilm phenotype marker |
| `frdB` | *Haemophilus influenzae* | B | Fumarate reductase — major literature gap |

---

## ⚙️ Pipeline Architecture

```
NCBI Entrez API
      │
      ▼
 fetch_and_align.py          ← Phase 0 entry (sequence acquisition + alignment)
      │
      ├─── MAFFT (local)     ← Multiple sequence alignment
      │
      ├─── Conservation Analysis
      │         └── Sliding window (20–25 nt)
      │             Criteria: conservation ≥ 85%, gap ≤ 20%
      │
      ▼
 gfet_probe_pipeline.py      ← Phase 0 full pipeline (v3, ViruScope integrated)
      │
      ├─── ViruScope Scoring (primer3-py)
      │         ├── Melting Temperature (NN method, SantaLucia 1998)
      │         ├── GC Content
      │         ├── Hairpin ΔG (no-fold criterion)
      │         └── Homodimer ΔG
      │
      └─── Outputs per target:
                ├── conserved_report.txt
                ├── <gene>_probes_scored.tsv
                ├── <gene>_viroscope_probes.fasta  ← ViruScope input
                └── aligned.fasta                  ← Geneious / AliView
```

---

## 📐 Probe Design Criteria

| Parameter | Threshold |
|-----------|-----------|
| Probe length | 20–25 nt |
| GC content | 40–60% |
| Melting temperature (Tm) | 55–72 °C |
| Hairpin ΔG | > −2.0 kcal/mol (no significant fold) |
| Homodimer ΔG | > −6.0 kcal/mol |
| Conservation (alignment) | ≥ 85% per position |
| Gap frequency | ≤ 20% |

Thermodynamic parameters use **Nearest-Neighbour method** (SantaLucia 1998) at 37 °C, 50 mM Na⁺, 250 nM oligo concentration — conditions representative of GFET biosensor operation.

---

## 🗂️ Repository Structure

```
.
├── fetch_and_align.py                   # Phase 0 — acquisition + alignment only
├── gfet_probe_pipeline.py               # Phase 0 — full pipeline with ViruScope scoring (v3)
├── tabela_targets_probes_updated.xlsx   # Targets and probe summary table
├── The_GFET_Computational_Blueprint.pdf # Full thesis document
├── Thesis_Plan_PG45861.pdf              # Thesis plan
├── Protocolo_Fase0_Alinhamento.docx     # Phase 0 lab protocol (PT)
└── documentacao_tecnica_fetch_and_align.docx  # Technical documentation (PT)
```

Output files (generated at runtime, not committed):
```
alignments/
└── <gene>/
    ├── raw_sequences.fasta
    ├── aligned.fasta
    ├── conserved_report.txt
    ├── <gene>_probes_scored.tsv
    └── <gene>_viroscope_probes.fasta
```

---

## 🚀 Getting Started

### Prerequisites

- Python ≥ 3.9
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (Windows: place at `MAFFT\mafft-7.526-win64-signed\mafft-win\mafft.bat`)
- Python packages:

```bash
pip install biopython primer3-py pandas
```

### Run the full pipeline

```bash
# All 6 targets
python gfet_probe_pipeline.py

# Specific target(s)
python gfet_probe_pipeline.py nuc lytA

# Skip NCBI fetch (use existing FASTA files)
python gfet_probe_pipeline.py --skip-fetch
```

### Output

A consolidated summary of the top 5 probes per target is saved to:
```
alignments/FINAL_PROBES_SUMMARY.tsv
```

---

## 🔁 Downstream Workflow (Phase 1)

After probe selection, the recommended validation pipeline is:

1. **Visual alignment review** — Geneious or AliView (`aligned.fasta`)
2. **ViruScope** — import `<gene>_viroscope_probes.fasta` for in-silico primer scoring
3. **NUPACK** ([nupack.org](https://www.nupack.org)) — secondary structure: ΔG < −12 kcal/mol, ensemble defect < 0.10
4. **BLAST** ([ncbi.nlm.nih.gov/blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) — cross-reactivity check (exclude target organism)
5. **3dDNA → AMBER → docking** — 3D modelling and molecular dynamics

---

## 📚 Key References

- Purwidyantri, A. et al. (2021). GFET biosensor for *S. aureus nuc* gene detection.
- Appaturi, J.N. et al. (2013). Electrochemical biosensor for *N. meningitidis rmpM*.
- SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS*, 95(4), 1460–1465.

---

## 📄 License

Academic use only. All rights reserved — University of Minho, 2025–2026.
