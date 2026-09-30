"""
Debye-screening analysis of the graphene + ssDNA probe models.

In an electrolyte-gated GFET only the charge lying within roughly one Debye length
(lambda_D) of the graphene surface contributes appreciably to the field-effect signal;
charge further away is screened by the mobile ions of the buffer. For each complex built
by `build_graphene_complex.py`, this script measures the height of every probe atom above
the graphene plane and reports the fraction of the probe that falls inside lambda_D.

The Debye lengths are the values reported by Purwidyantri et al., Biosensors 11(1):24
(2021), who calculated them with the Debye-Huckel approximation for the buffers used in
their DNA-hybridisation experiments on graphene transistors:

    1x   PBS   ->  0.76 nm  ( 7.6 A)   osmolarity and ion content close to body fluids
    0.1x PBS   ->  2.41 nm  (24.1 A)
    0.01x PBS  ->  7.61 nm  (76.1 A)

Outputs
    results/debye_coverage.csv   one row per probe: height above the sheet and covered fraction
    results/debye_by_class.csv   the same aggregated by Boltz-2 confidence class

Usage:
    python scripts/debye_analysis.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
COMPLEX_DIR = BASE_DIR / "results" / "structures" / "complexes"
SUMMARY = BASE_DIR / "results" / "complexes_summary.csv"

# Debye lengths, in Angstrom (Purwidyantri et al., Biosensors 11(1):24, 2021)
BUFFERS = {"1x PBS": 7.6, "0.1x PBS": 24.1, "0.01x PBS": 76.1}


def atom_heights(pdb_path):
    """Heights (A) of every chain-D (probe) atom above the mean chain-G (graphene) plane."""
    z_gra, z_dna = [], []
    with open(pdb_path) as fh:
        for line in fh:
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            z = float(line[46:54])
            (z_gra if line[21] == "G" else z_dna).append(z)
    return np.array(z_dna) - np.mean(z_gra)


def main():
    meta = pd.read_csv(SUMMARY).set_index("probe")

    rows = []
    for pdb in sorted(COMPLEX_DIR.glob("graphene_*.pdb")):
        probe = pdb.stem.split("_")[1]
        h = atom_heights(pdb)
        row = {"probe": probe,
               "gene": meta.at[probe, "gene"],
               "atoms": len(h),
               "height_above_sheet": round(float(h.max()), 2),
               "confidence": meta.at[probe, "confidence"],
               "plddt": meta.at[probe, "plddt"],
               "boltz_class": meta.at[probe, "boltz_class"]}
        for name, lam in BUFFERS.items():
            row[f"f_{name}"] = round(100.0 * float((h <= lam).mean()), 1)
        rows.append(row)

    cov = pd.DataFrame(rows).sort_values(["gene", "probe"])
    # ";" so that the tables open in columns in Excel with a Portuguese locale
    cov.to_csv(BASE_DIR / "results" / "debye_coverage.csv", index=False, sep=";")

    # aggregate by Boltz-2 confidence class
    order = ["GOOD", "MODERATE", "LOW"]
    agg = []
    for cls in order:
        sub = cov[cov["boltz_class"] == cls]
        if sub.empty:
            continue
        rec = {"boltz_class": cls, "n": len(sub),
               "height_median": round(sub["height_above_sheet"].median(), 1),
               "height_min": round(sub["height_above_sheet"].min(), 1),
               "height_max": round(sub["height_above_sheet"].max(), 1)}
        for name in BUFFERS:
            col = sub[f"f_{name}"]
            rec[f"{name} median %"] = round(col.median(), 1)
            rec[f"{name} min %"] = round(col.min(), 1)
            rec[f"{name} max %"] = round(col.max(), 1)
            rec[f"{name} n>=50%"] = int((col >= 50).sum())
        agg.append(rec)
    agg = pd.DataFrame(agg)
    agg.to_csv(BASE_DIR / "results" / "debye_by_class.csv", index=False, sep=";")

    print(cov.to_string(index=False))
    print()
    print(agg.to_string(index=False))
    print()
    for name in BUFFERS:
        r = cov[f"f_{name}"].corr(cov["confidence"])
        print(f"  Pearson r(coverage, Boltz confidence) at {name:9s}: {r:+.2f}")


if __name__ == "__main__":
    main()
