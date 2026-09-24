"""
Build starting-pose graphene + ssDNA probe complexes.

Places a Boltz-2-predicted probe structure (CIF) above a GOPY pristine-graphene sheet
(PDB): the probe is translated so that its centroid lies over the centre of the sheet
in x-y, and so that its lowest atom sits `gap` Angstrom above the graphene plane in z.
The probe's internal conformation and orientation are left untouched.

The default separation, 3.2 A, is the equilibrium physisorption height of DNA nucleobases
on graphene from van der Waals-corrected DFT (3.16-3.22 A; Tayo, Walkup & Caliskan, AIP
Advances 13:085213, 2023; see also Lee et al., J. Phys. Chem. C 117:13435, 2013), so the
model starts at the literature contact distance rather than an arbitrary one.

This is an initial physical placement (a starting geometry for subsequent docking or
force-field refinement), not a docked or energy-minimised pose.

The sheet must be large enough to host the probe: the script reports the probe
footprint and refuses to write a complex in which probe atoms would overhang the sheet
edges (use a larger sheet, or --allow-overhang to override).

Usage:
    # one model
    python scripts/build_graphene_complex.py <probe.cif> <graphene.pdb> <out.pdb>
                                             [--gap 3.2] [--allow-overhang]

    # every probe of a directory, plus a summary table
    python scripts/build_graphene_complex.py <probes_dir/> <graphene.pdb> <out_dir/>
                                             [--meta docs/boltz2_shortlist_ranked.csv]
"""
import argparse
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

BASE_DIR = Path(__file__).resolve().parent.parent
DEFAULT_META = BASE_DIR / "docs" / "boltz2_shortlist_ranked.csv"


def build_one(probe_path, graphene_path, out_path, gap, allow_overhang, quiet=False):
    """Place one probe on the sheet, write the merged PDB, return its measurements."""
    gra = PDBParser(QUIET=True).get_structure("gra", str(graphene_path))
    prb = MMCIFParser(QUIET=True).get_structure("prb", str(probe_path))

    g = np.array([a.coord for a in gra.get_atoms()])
    prb_atoms = list(prb.get_atoms())
    p = np.array([a.coord for a in prb_atoms])

    # centre the probe over the sheet (x, y) and lift it `gap` above the plane (z)
    shift = np.array([g[:, 0].mean() - p[:, 0].mean(),
                      g[:, 1].mean() - p[:, 1].mean(),
                      (g[:, 2].mean() + gap) - p[:, 2].min()])
    for a in prb_atoms:
        a.set_coord(a.coord + shift)
    p = np.array([a.coord for a in prb_atoms])

    # the probe must sit within the sheet footprint
    outside = int(((p[:, 0] < g[:, 0].min()) | (p[:, 0] > g[:, 0].max()) |
                   (p[:, 1] < g[:, 1].min()) | (p[:, 1] > g[:, 1].max())).sum())
    dmin = np.linalg.norm(p[:, None, :] - g[None, :, :], axis=-1).min()  # closest contact

    if not quiet:
        print(f"  sheet    : {len(g)} atoms, {np.ptp(g[:,0]):.1f} x {np.ptp(g[:,1]):.1f} A")
        print(f"  probe    : {len(p)} atoms, footprint {np.ptp(p[:,0]):.1f} x {np.ptp(p[:,1]):.1f} A")
        print(f"  contact  : closest probe-sheet atom pair {dmin:.2f} A")
        print(f"  overhang : {outside} probe atoms outside the sheet")
    if outside and not allow_overhang:
        raise SystemExit("  ! probe overhangs the sheet - use a larger sheet (or --allow-overhang)")

    # merge into a single model: graphene as chain G, probe as chain D
    model = gra[0]
    list(model.get_chains())[0].id = "G"
    chain = list(prb[0].get_chains())[0]
    chain.detach_parent()
    chain.id = "D"
    model.add(chain)

    io = PDBIO()
    io.set_structure(gra)
    io.save(str(out_path))
    if not quiet:
        print(f"  -> {out_path}  (gap {gap} A, chains G=graphene, D=probe)")

    return {"probe_atoms": len(p), "total_atoms": len(g) + len(p),
            "footprint_x": round(float(np.ptp(p[:, 0])), 1),
            "footprint_y": round(float(np.ptp(p[:, 1])), 1),
            "height": round(float(np.ptp(p[:, 2])), 1),
            "min_contact": round(float(dmin), 2),
            "overhang": outside}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("probe", help="probe CIF, or a directory of probe CIFs (batch mode)")
    ap.add_argument("graphene")
    ap.add_argument("out", help="output PDB, or output directory in batch mode")
    ap.add_argument("--gap", type=float, default=3.2,
                    help="vertical separation probe-to-sheet, in Angstrom (default 3.2, the "
                         "equilibrium nucleobase-graphene physisorption height from vdW-DFT)")
    ap.add_argument("--allow-overhang", action="store_true",
                    help="write the complex even if the probe extends beyond the sheet")
    ap.add_argument("--meta", default=str(DEFAULT_META),
                    help="CSV with the Boltz-2 metrics, joined into the batch summary")
    args = ap.parse_args()

    probe = Path(args.probe)
    if not probe.is_dir():
        build_one(probe, args.graphene, args.out, args.gap, args.allow_overhang)
        return

    # batch mode: every CIF of the directory, plus a summary table
    import pandas as pd

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for cif in sorted(probe.glob("*.cif")):
        # file names follow <probe>_<organism>_<gene>.cif
        parts = cif.stem.split("_")
        name, gene = parts[0], parts[-1]
        out_path = out_dir / f"graphene_{name}_{gene}.pdb"
        rec = build_one(cif, args.graphene, out_path, args.gap, args.allow_overhang, quiet=True)
        rows.append({"probe": name, "gene": gene, **rec})
        print(f"  {name:6s} {gene:5s} footprint {rec['footprint_x']:5.1f} x "
              f"{rec['footprint_y']:5.1f} A, contact {rec['min_contact']:.2f} A -> {out_path.name}")

    summary = pd.DataFrame(rows)
    meta_path = Path(args.meta)
    if meta_path.exists():
        meta = pd.read_csv(meta_path)
        meta["probe"] = meta["probe_id"].astype(str).str.split("_").str[0]
        keep = [c for c in ["probe", "confidence", "plddt", "quality"] if c in meta.columns]
        summary = summary.merge(meta[keep], on="probe", how="left")
        summary = summary.rename(columns={"quality": "boltz_class"})

    csv_path = out_dir.parent.parent / "complexes_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"\n  {len(summary)} models built, summary -> {csv_path}")


if __name__ == "__main__":
    main()
