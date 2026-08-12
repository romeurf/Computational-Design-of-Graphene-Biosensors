# Results — 3D structures (probe–graphene modelling)

3D structures for the ssDNA probe (aptamer)–graphene modelling stage.

## Contents

| Path | What it is |
|---|---|
| `structures/graphene_PG_12nm.pdb` | Pristine graphene sheet, 121.6 × 119.8 Å (5700 C atoms) |
| `structures/aptamers/*.cif` | 30 Boltz-2-predicted ssDNA probe structures (5 per gene) |
| `structures/complexes/*.pdb` | Graphene + probe models (chain `G` = graphene, `D` = probe) |
| `boltz2_all_results.zip` | Raw Boltz-2 output (CIF + confidence JSON + PAE NPZ) |

## How each piece was made

**1. Graphene sheet** — [GOPY](https://github.com/Iourarum/GOPY) (Muraru, Burns & Ioniță,
*SoftwareX* 12:100586, 2020):

```
python GOPY.py generate_PG 120 120 graphene_PG_12nm.pdb
```

The size is not arbitrary: the largest in-plane extent of the predicted probes reaches 82.9 Å, so
a 4 × 4 nm sheet would be too small for 16 of the 30 probes. At 121.6 × 119.8 Å all 30 fit, the
largest with a 19.3 Å margin per side.

*(Note: this is the graphene-modelling GOPY, not `go-python/gopy`, an unrelated Go↔Python
bindings generator.)*

**2. Probe structures** — Boltz-2, via `colab_boltz2_batch.ipynb`; the `_model_0.cif` of each
prediction was extracted from the raw archive.

**3. Complexes** — [`scripts/build_graphene_complex.py`](../scripts/build_graphene_complex.py):

```
python scripts/build_graphene_complex.py results/structures/aptamers/p5972_Paer_algD.cif \
       results/structures/graphene_PG_12nm.pdb \
       results/structures/complexes/graphene_p5972_algD.pdb --gap 5.0
```

The probe is translated (rigid body, no distortion) so that it is centred over the sheet in x–y
and its lowest atom sits 5 Å above the graphene plane; the two are merged into one PDB as chains
`G` and `D`. The script refuses to write a model in which probe atoms would overhang the sheet.

On the 5 Å separation: vdW-corrected DFT places physisorbed nucleobases at **3.16–3.22 Å** above
graphene (Tayo, Walkup & Caliskan, *AIP Advances* 13:085213, 2023; see also Lee et al.,
*J. Phys. Chem. C* 117:13435, 2013). The starting pose is set deliberately a little above that
equilibrium contact distance, so the geometry is clash-free and the final separation is left to
docking/refinement instead of being fixed by construction. It is the `--gap` parameter.

> These are **starting geometries** — an initial physical placement for subsequent docking or
> force-field refinement. They are not docked, energy-minimised or equilibrated structures, and no
> interaction energy is computed.

Built so far, for the two highest-confidence probes: `graphene_p5972_algD.pdb` (6070 atoms) and
`graphene_p7716_frdB.pdb` (6229 atoms).
