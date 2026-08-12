"""
Build a starting-pose graphene + ssDNA probe complex.

Places a Boltz-2-predicted probe structure (CIF) above a GOPY pristine-graphene sheet
(PDB): the probe is translated so that its centroid lies over the centre of the sheet
in x-y, and so that its lowest atom sits `gap` Angstrom above the graphene plane in z.
The probe's internal conformation and orientation are left untouched.

This is an initial physical placement (a starting geometry for subsequent docking or
force-field refinement), not a docked or energy-minimised pose.

The sheet must be large enough to host the probe: the script reports the probe
footprint and refuses to write a complex in which probe atoms would overhang the sheet
edges (use a larger sheet, or --allow-overhang to override).

Usage:
    python scripts/build_graphene_complex.py <probe.cif> <graphene.pdb> <out.pdb>
                                             [--gap 5.0] [--allow-overhang]
"""
import argparse
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

ap = argparse.ArgumentParser()
ap.add_argument("probe"); ap.add_argument("graphene"); ap.add_argument("out")
ap.add_argument("--gap", type=float, default=5.0,
                help="vertical separation probe-to-sheet, in Angstrom (default 5.0)")
ap.add_argument("--allow-overhang", action="store_true",
                help="write the complex even if the probe extends beyond the sheet")
args = ap.parse_args()

gra = PDBParser(QUIET=True).get_structure("gra", args.graphene)
prb = MMCIFParser(QUIET=True).get_structure("prb", args.probe)

g = np.array([a.coord for a in gra.get_atoms()])
prb_atoms = list(prb.get_atoms())
p = np.array([a.coord for a in prb_atoms])

# centre the probe over the sheet (x, y) and lift it `gap` above the plane (z)
shift = np.array([g[:, 0].mean() - p[:, 0].mean(),
                  g[:, 1].mean() - p[:, 1].mean(),
                  (g[:, 2].mean() + args.gap) - p[:, 2].min()])
for a in prb_atoms:
    a.set_coord(a.coord + shift)
p = np.array([a.coord for a in prb_atoms])

# the probe must sit within the sheet footprint
outside = int(((p[:, 0] < g[:, 0].min()) | (p[:, 0] > g[:, 0].max()) |
               (p[:, 1] < g[:, 1].min()) | (p[:, 1] > g[:, 1].max())).sum())
print(f"  sheet    : {len(g)} atoms, {np.ptp(g[:,0]):.1f} x {np.ptp(g[:,1]):.1f} A")
print(f"  probe    : {len(p)} atoms, footprint {np.ptp(p[:,0]):.1f} x {np.ptp(p[:,1]):.1f} A")
print(f"  overhang : {outside} probe atoms outside the sheet")
if outside and not args.allow_overhang:
    raise SystemExit("  ! probe overhangs the sheet - use a larger sheet (or --allow-overhang)")

# merge into a single model: graphene as chain G, probe as chain D
model = gra[0]
list(model.get_chains())[0].id = "G"
chain = list(prb[0].get_chains())[0]
chain.detach_parent()
chain.id = "D"
model.add(chain)

io = PDBIO(); io.set_structure(gra); io.save(args.out)
print(f"  -> {args.out}  (gap {args.gap} A, chains G=graphene, D=probe)")
