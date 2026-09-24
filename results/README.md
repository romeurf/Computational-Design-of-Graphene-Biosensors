# Results — 3D structures (probe–graphene modelling)

3D structures for the ssDNA probe (aptamer)–graphene modelling stage, and the measurements taken
from them.

## Contents

| Path | What it is |
|---|---|
| `structures/graphene_PG_12nm.pdb` | Pristine graphene sheet, 121.6 × 119.8 Å (5700 C atoms) |
| `structures/aptamers/*.cif` | 30 Boltz-2-predicted ssDNA probe structures (5 per gene) |
| `structures/complexes/*.pdb` | 30 graphene + probe models (chain `G` = graphene, `D` = probe) |
| `boltz2_all_results.zip` | Raw Boltz-2 output (CIF + confidence JSON + PAE NPZ) |
| `complexes_summary.csv` | Per model: atom counts, footprint, height, closest contact, Boltz-2 confidence |
| `debye_coverage.csv` | Per model: fraction of the probe inside the Debye length, per PBS dilution |
| `debye_by_class.csv` | The same, aggregated by Boltz-2 confidence class |
| `complex_p5972_algD.png`, `complex_p7716_frdB.png` | Side + top renders of two representative models |

## How each piece was made

**1. Graphene sheet** — [GOPY](https://github.com/Iourarum/GOPY) (Muraru, Burns & Ioniță,
*SoftwareX* 12:100586, 2020). GOPY is a standalone GPL-licensed Python script, not a package on
PyPI, so it is used by cloning the author's repository and running the script directly:

```
python GOPY.py generate_PG 120 120 graphene_PG_12nm.pdb
```

`generate_PG` builds a pristine graphene sheet (PG) of the requested size in Ångström. The size is
not arbitrary: the largest in-plane extent among the 30 predicted probes is 82.9 Å, so a 40 × 40 Å
sheet would be too small for 16 of them. At 121.6 × 119.8 Å all 30 fit inside the sheet, the
largest with a 19.3 Å margin per side.

*(Note: this is the graphene-modelling GOPY, not `go-python/gopy` — an unrelated Go↔Python bindings
generator — nor the `gopy` package on PyPI, which is a data-structures teaching library.)*

**2. Probe structures** — Boltz-2, via `colab_boltz2_batch.ipynb`; the `_model_0.cif` of each
prediction was extracted from the raw archive.

**3. Complexes** — [`scripts/build_graphene_complex.py`](../scripts/build_graphene_complex.py):

```
# one model
python scripts/build_graphene_complex.py results/structures/aptamers/p5972_Paer_algD.cif \
       results/structures/graphene_PG_12nm.pdb \
       results/structures/complexes/graphene_p5972_algD.pdb

# the whole panel in one call, plus complexes_summary.csv
python scripts/build_graphene_complex.py results/structures/aptamers \
       results/structures/graphene_PG_12nm.pdb results/structures/complexes
```

The probe is translated (rigid body, no distortion) so that it is centred over the sheet in x–y and
its lowest atom sits `--gap` Å above the graphene plane; the two are merged into one PDB as chains
`G` (graphene) and `D` (DNA). The script reports the probe footprint and refuses to write a model
in which probe atoms would overhang the sheet.

The default `--gap` is **3.2 Å**, the equilibrium physisorption height of DNA nucleobases on
graphene from van der Waals-corrected DFT (3.16–3.22 Å; Tayo, Walkup & Caliskan, *AIP Advances*
13:085213, 2023; see also Lee et al., *J. Phys. Chem. C* 117:13435, 2013). The parameter is set
from that literature value rather than chosen arbitrarily. Across the 30 models the closest
probe–sheet atom pair falls between 3.20 and 3.47 Å (median 3.25 Å), with no steric clash and no
overhang.

> These are **starting geometries** — an initial physical placement for subsequent docking or
> force-field refinement. They are not docked, energy-minimised or equilibrated structures, and no
> interaction energy is computed.

**4. Debye-screening measurements** — [`scripts/debye_analysis.py`](../scripts/debye_analysis.py)
measures the height of every probe atom above the graphene plane and reports how much of each probe
lies within the Debye length λ_D of the three PBS dilutions measured on graphene FETs by
Purwidyantri et al., *Biosensors* 11:120 (2021): 0.76 nm (1×), 2.41 nm (0.1×) and 7.61 nm (0.01×).

**5. Renders** — [`scripts/render_complex.py`](../scripts/render_complex.py) produces the side and
top views, with the probe–sheet separation annotated.

## What the measurements show

| Boltz-2 class | *n* | Median height (Å) | Within λ_D at 1× | at 0.1× | at 0.01× |
|---|---|---|---|---|---|
| GOOD | 4 | 41.4 | 5.6% | 51.8% | 100% |
| MODERATE | 11 | 44.3 | 6.7% | 50.9% | 100% |
| LOW | 15 | 33.7 | 7.5% | 76.9% | 100% |

Coverage is essentially independent of the 3D confidence (Pearson *r* = −0.19, −0.26 and −0.22 at
1×, 0.1× and 0.01× PBS), so interface geometry is a selection axis of its own. At 1× PBS — the
isotonic condition — no probe exceeds 15% coverage; at 0.01× PBS 29 of the 30 models are entirely
inside λ_D and the criterion no longer distinguishes candidates.
