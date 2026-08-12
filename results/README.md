# Results — 3D structures (probe–graphene modelling)

3D structures for the ssDNA probe (aptamer)–graphene modelling stage.

## Contents

| Path | What it is |
|---|---|
| `structures/graphene_PG_4nm.pdb` | Pristine graphene sheet, 4 × 4 nm (646 C atoms) |
| `structures/aptamers/*.cif` | 30 Boltz-2-predicted ssDNA probe structures (5 per gene) |
| `structures/complexes/*.pdb` | Graphene + probe starting-pose models |
| `boltz2_all_results.zip` | Raw Boltz-2 output (CIF + confidence JSON + PAE NPZ) |

## How each piece was made

- **Graphene** — [GOPY](https://github.com/Iourarum/GOPY) (Muraru & Ionita, *SoftwareX* 2020):
  ```
  python GOPY.py generate_PG 40 40 graphene_PG.pdb
  ```
  *(This is the graphene-modelling GOPY, not `go-python/gopy`, which is an unrelated
  Go↔Python bindings generator.)*
- **Probe structures** — Boltz-2, via `colab_boltz2_batch.ipynb`.
- **Complexes** — `scripts/build_graphene_complex.py`, which centres a probe over the graphene
  sheet and lifts it 5 Å above the plane:
  ```
  python scripts/build_graphene_complex.py <probe.cif> structures/graphene_PG_4nm.pdb <out.pdb>
  ```
  This is an **initial physical placement** (a starting pose for docking/refinement), not a
  docked or energy-minimised structure.

The two complexes built so far use the two highest-confidence probes: `p5972` (*algD*) and
`p7716` (*frdB*).
