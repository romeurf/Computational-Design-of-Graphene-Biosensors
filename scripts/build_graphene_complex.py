"""
Build a starting-pose graphene + ssDNA (aptamer/probe) complex.

Places a Boltz-2-predicted probe structure (CIF) above a GOPY pristine-graphene
sheet (PDB): the probe is centred over the sheet and lifted so its lowest atom sits
`gap` Angstrom above the graphene plane. This is an initial physical placement for
subsequent docking/refinement, not a docked pose.

Usage:
    python scripts/build_graphene_complex.py <aptamer.cif> <graphene.pdb> <out.pdb> [gap_A]
"""
import sys
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

apt_path, gra_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
gap = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0

gra = PDBParser(QUIET=True).get_structure("gra", gra_path)
apt = MMCIFParser(QUIET=True).get_structure("apt", apt_path)

g = np.array([a.coord for a in gra.get_atoms()])
apt_atoms = list(apt.get_atoms())
a = np.array([at.coord for at in apt_atoms])

# centre the probe over the sheet (x, y) and lift it `gap` A above the plane (z)
tran = np.array([g[:, 0].mean() - a[:, 0].mean(),
                 g[:, 1].mean() - a[:, 1].mean(),
                 (g[:, 2].mean() + gap) - a[:, 2].min()])
for at in apt_atoms:
    at.set_coord(at.coord + tran)

# merge into one model: graphene chain "G", probe chain "D"
gmodel = gra[0]
list(gmodel.get_chains())[0].id = "G"
apt_chain = list(apt[0].get_chains())[0]
apt_chain.detach_parent()
apt_chain.id = "D"
gmodel.add(apt_chain)

io = PDBIO()
io.set_structure(gra)
io.save(out_path)
print(f"wrote {out_path}: graphene ({g.shape[0]} atoms) + probe ({a.shape[0]} atoms), gap {gap} A")
