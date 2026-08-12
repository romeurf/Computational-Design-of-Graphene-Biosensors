# Results — 3D structures (aptamer–graphene modelling)

3D structures for the aptamer/probe–graphene modelling stage.

## Graphene
- `structures/graphene_PG_4nm.pdb` — pristine graphene sheet (4 × 4 nm, 646 C atoms),
  generated with **GOPY** ([Iourarum/GOPY](https://github.com/Iourarum/GOPY), *SoftwareX* 2020):
  ```
  python GOPY.py generate_PG 40 40 graphene_PG.pdb
  ```
  (Note: this is the graphene-modelling GOPY, not `go-python/gopy`, which is an unrelated
  Go↔Python bindings generator.)

## Aptamer / probe
- 3D structures of the selected ssDNA probes, predicted with Boltz-2 (to be added from the
  Colab run — the CIF/PDB files).
