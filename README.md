# Computational Design of Graphene Biosensors

![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python&logoColor=white)
![University](https://img.shields.io/badge/University-Minho-darkgreen)

Computational pipeline for the *in silico* design of single-stranded DNA (ssDNA) capture probes
for **Graphene Field-Effect Transistor (GFET) biosensors**. Given a target gene, the pipeline goes
from public NCBI sequences to a small set of ranked, biophysically screened probe candidates and
prepares them for 3D structure validation.

Everything runs from a single script, [`pipeline.py`](pipeline.py).

---

## Documents

- 📄 **[Thesis (PDF)](Thesis.pdf)** — dissertation document (view / download)
- 📝 **[Thesis (DOCX)](Thesis.docx)** — editable version

*(The LaTeX source of the dissertation is under [`thesis/`](thesis/).)*

---

## Overview

```
NCBI Entrez  →  length clustering  →  MAFFT alignment  →  PPI conservation
     →  sliding-window candidates  →  primer3 thermodynamics  →  seqfold + No-fold
     →  organism-aware thresholds  →  quality ranking  →  Boltz-2 3D export
```

All thresholds are resolved through a layered scheme (**gene → species → organism type → global
default**), so AT-rich and GC-rich targets are screened by criteria appropriate to their own genome.

---

## Repository structure

```
.
├── pipeline.py                 Main entry point (all stages + CLI)
├── scripts/
│   ├── analysis.py             Exploratory analysis (diversity, rarefaction, seqfold threshold)
│   ├── boltz2_predict.py       Standalone local Boltz-2 predictor
│   └── nucleofold_predict.py   NucleoFold3D wrapper (WSL/Linux)
├── colab_boltz2_batch.ipynb    Boltz-2 batch 3D prediction on Colab (GPU)
├── data/
│   ├── sequencias_iplex.xlsx   Reference probe panel (job_name, sequence, species)
│   └── species_params.yaml     User-defined species profiles (created at runtime)
├── docs/                       Reference tables, deliverables, analysis outputs, notes
├── thesis/                     LaTeX dissertation (UMinho template, XeLaTeX)
└── output/                     Generated at runtime (git-ignored)
```

> **MAFFT** is used as an external aligner. Install it and make sure `mafft` is on the `PATH`
> (the pipeline also auto-detects a local `MAFFT/` folder in the repository root if one exists).

---

## Requirements

```bash
pip install biopython primer3-py pandas pyyaml seqfold openpyxl
# for scripts/analysis.py:
pip install numpy matplotlib scipy scikit-learn
```

External tools:
- **MAFFT** — multiple sequence alignment (`mafft` on `PATH`).
- **Boltz-2** — 3D structure prediction, run on a GPU (via the Colab notebook, or `pip install boltz` locally).

---

## Pipeline stages

Each stage is a function in `pipeline.py`:

| Stage | Function(s) | What it does |
|---|---|---|
| 1. Retrieval | `fetch_sequences` | NCBI Entrez `esearch`/`efetch` (batched); primary → alternative → RefSeq fallback queries |
| 2. Length clustering | `select_length_cluster` | Keeps the dominant length band (±25%), dropping partial/oversized records |
| 3. Alignment | `align_mafft` | MAFFT `--auto` (single-sequence mode when only one record) |
| 4. Conservation | `_pairwise_identity`, `candidate_windows` | Per-column **PPI** (Percentage of Pairwise Identity, IUPAC-aware; ported from ViruScope) |
| 5. Candidates | `candidate_windows` | Sliding window (18–28 nt, step 3), gap and conservation filters, consensus sequence |
| 6. Thermodynamics | `score_probe` | primer3-py: Tm, GC, hairpin ΔG, homodimer ΔG ([Na⁺]=50 mM, [oligo]=250 nM, 37 °C) |
| 7. Structure | `run_seqfold_probe`, `nofold_score` | seqfold MFE at 37 °C; continuous **No-fold** score (logistic map of the folding energies) |
| 8. Parameters | `cfg`, `cfg_species`, `_infer_type` | Layered threshold resolution + interactive profiles for new species |
| 9. Selection | `probe_quality`, `export_colab_inputs` | Rank by quality = ½·(PPI + No-fold/100); export top-N/gene |
| 10. 3D export | `export_colab_from_csv`, `merge_boltz_results` | Boltz-2 YAML/zip input; merge results into a ranked shortlist |

### Configuration (in `pipeline.py`)

- `DEFAULTS` — global thresholds (length 18–28 nt, Tm 53–72 °C, GC 40–60%, hairpin ≥ −2.0,
  homodimer ≥ −5.0, PPI ≥ 0.85, gap ≤ 20%, seqfold MFE ≥ −2.6 kcal/mol).
- `TYPE_DEFAULTS` — five organism profiles (bacteria / virus / fungus / protozoa / host).
- `SPECIES_PARAMS` — per-species overrides tuned by genome GC content.
- `TARGETS` — the six target genes and their NCBI queries.
- `Probe` — dataclass holding every per-probe metric.

---

## Usage

```bash
# Full run over all targets (100 NCBI sequences per gene), export top-5/gene for Boltz-2
python pipeline.py --max-seqs 100 --colab 5

# Include the reference (IPLEX) panel, scored with the same metrics, in one run
python pipeline.py --max-seqs 100 --with-reference --colab 5

# Merge Boltz-2 results (downloaded from Colab) into a ranked shortlist
python pipeline.py --merge-boltz output/colab_boltz2/boltz2_results_summary.csv

# Export a curated set for external use (top-N per gene by quality)
python pipeline.py --export-docking 25

# Define a profile for a new species interactively (type + thresholds, with suggestions)
python pipeline.py --define-species "Influenza A H3N2"
```

Main flags:

| Flag | Effect |
|---|---|
| `--max-seqs N` | NCBI sequences to fetch per gene (default 20) |
| `--colab N` | Full run, then export the top-N/gene (by quality) as Boltz-2 input |
| `--with-reference` | Also score the reference panel (`data/sequencias_iplex.xlsx`) |
| `--export-colab N` / `--export-colab-iplex N` | Re-export Boltz-2 input from an existing `FINAL_PROBES_ALL.csv` |
| `--merge-boltz CSV` | Merge Boltz-2 results into `boltz2_shortlist_ranked.csv` |
| `--export-docking [N]` | Export `job_name,sequence,species` (+features); with `N`, top-N/gene |
| `--define-species NAME` | Interactive per-species profile → `data/species_params.yaml` |
| `--assay-temp C` | Derive the Tm window from the hybridisation assay temperature |
| `--no-cluster` / `--cluster-tol T` | Control the length-clustering step |
| `--no-seqfold` | Skip the secondary-structure stage |

### Outputs (in `output/`)

| File | Description |
|---|---|
| `FINAL_PROBES_ALL.csv` | All probes (own + optional reference) with every metric |
| `alignments/<gene>/` | MAFFT alignment, per-gene scored TSV and FASTA |
| `colab_boltz2/boltz2_inputs.zip` | Boltz-2 input (one YAML per probe) for the Colab notebook |
| `colab_boltz2/boltz2_shortlist_ranked.csv` | Ranked shortlist after merging 3D results |

---

## 3D structure validation

The selected probes are exported as Boltz-2 inputs and predicted on a GPU:

1. Run `pipeline.py ... --colab N` to produce `output/colab_boltz2/boltz2_inputs.zip`.
2. Open [`colab_boltz2_batch.ipynb`](colab_boltz2_batch.ipynb) in Google Colab (GPU runtime), upload the zip, run it.
3. Download the results and run `pipeline.py --merge-boltz <results.csv>`.

Standalone local predictors are also available under `scripts/` (`boltz2_predict.py`,
`nucleofold_predict.py`).

---

## Exploratory analysis

```bash
python scripts/analysis.py
```

Writes figures and tables to `output/figures/` and `output/analysis/`: seqfold ΔG distribution and
threshold sweep, per-gene sizes and normality, alignment-free k-mer diversity, and a rarefaction
analysis used to justify the sampling depth.

---

## Thesis

The accompanying MSc dissertation (UMinho template, compiled with **XeLaTeX**) lives in
[`thesis/`](thesis/).

---

*Academic use — University of Minho.*
