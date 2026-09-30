"""
Debye-screening analysis of the graphene + ssDNA probe models.

In an electrolyte-gated GFET only the charge lying within roughly one Debye length
(lambda_D) of the graphene surface contributes appreciably to the field-effect signal;
charge further away is screened by the mobile ions of the buffer. DNA carries its charge on
the phosphates (one negative charge each), so besides the fraction of atoms inside lambda_D
the script counts the phosphorus atoms inside it.

Three geometries are measured:
  1. rigid placement - the complexes built by `build_graphene_complex.py` (the Boltz-2
     orientation, lowest atom at the contact distance);
  2. anchored probe - the probe tethered by its 3' end, as when it is immobilised through a
     3'-amino linker and PBSE (Campos et al. 2019; Purwidyantri et al. 2021). The probe takes
     N_ORIENTATIONS random orientations; in each, its lowest atom is placed at the contact
     distance, and the orientation is kept only if that atom belongs to the 3'-terminal
     nucleotide, i.e. if the probe touches the sheet with its 3' end. The median and maximum
     over the kept orientations are reported. No linker length is added, so the 3' end sits
     as close to the sheet as it can;
  3. anchored duplex - the same for the probe hybridised to its complementary target (Boltz-2
     duplex: chain A = probe, chain B = target), counting the target's phosphates, i.e. the
     charge that hybridisation brings inside lambda_D.

For heights above a flat sheet only the direction that ends up pointing away from the sheet
matters (a rotation about it changes no height), so an orientation is a direction drawn
uniformly on the sphere.

The Debye lengths are the values reported by Purwidyantri et al., Biosensors 11(1):24
(2021), who calculated them with the Debye-Huckel approximation for the buffers used in
their DNA-hybridisation experiments on graphene transistors:

    1x   PBS   ->  0.76 nm  ( 7.6 A)   osmolarity and ion content close to body fluids
    0.1x PBS   ->  2.41 nm  (24.1 A)
    0.01x PBS  ->  7.61 nm  (76.1 A)

Final selection: in each gene, the probes on the Pareto front of Boltz-2 confidence and
anchored phosphates inside lambda_D at 1x PBS (no other probe of the gene is at least as
good on both and better on one).

Outputs
    results/debye_coverage.csv   one row per probe
    results/debye_by_class.csv   the same aggregated by Boltz-2 confidence class

Usage:
    python scripts/debye_analysis.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS = BASE_DIR / "results"
COMPLEX_DIR = RESULTS / "structures" / "complexes"
PROBE_DIR = RESULTS / "structures" / "aptamers"
DUPLEX_DIR = RESULTS / "structures" / "duplexes"
SUMMARY = RESULTS / "complexes_summary.csv"

# Debye lengths, in Angstrom (Purwidyantri et al., Biosensors 11(1):24, 2021)
BUFFERS = {"1x PBS": 7.6, "0.1x PBS": 24.1, "0.01x PBS": 76.1}
GAP = 3.2              # contact distance, A (Tayo et al. 2023) - as in build_graphene_complex.py
N_ORIENTATIONS = 20000
SEED = 0


def complex_heights(pdb_path):
    """Heights (A) of every chain-D (probe) atom above the mean chain-G (graphene) plane,
    and which of those atoms are phosphorus."""
    z_gra, z_dna, is_p = [], [], []
    with open(pdb_path) as fh:
        for line in fh:
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            z = float(line[46:54])
            if line[21] == "G":
                z_gra.append(z)
            else:
                z_dna.append(z)
                is_p.append(line[76:78].strip() == "P")
    return np.array(z_dna) - np.mean(z_gra), np.array(is_p)


def up_directions(n, seed):
    """n directions drawn uniformly on the sphere (normalised Gaussian vectors)."""
    u = np.random.default_rng(seed).normal(size=(n, 3))
    return u / np.linalg.norm(u, axis=1, keepdims=True)


def load_model(cif_path):
    """Coordinates, phosphorus mask and chain of every atom, the mask of the 3'-terminal
    nucleotide of chain A (the probe) and the structure itself."""
    s = MMCIFParser(QUIET=True).get_structure("m", str(cif_path))
    atoms = list(s[0].get_atoms())
    xyz = np.array([a.coord for a in atoms], dtype=np.float32)
    is_p = np.array([a.element == "P" for a in atoms])
    chain = np.array([a.get_parent().get_parent().id for a in atoms])
    last = list(s[0]["A"])[-1]
    terminal = np.array([a.get_parent() is last for a in atoms])
    return xyz, is_p, chain, terminal, s


def anchored_heights(xyz, terminal, u):
    """Heights above the sheet (atoms x kept orientations), lowest atom at GAP; an
    orientation is kept when its lowest atom belongs to the 3'-terminal nucleotide."""
    h = xyz @ u.T.astype(np.float32)
    keep = terminal[h.argmin(axis=0)]
    h = h[:, keep]
    return h - h.min(axis=0) + GAP


def rise_per_bp(structure):
    """Mean rise per base pair (A) of a duplex: distance between the C1' midpoints of the
    first and last base pairs over the number of steps (B-DNA: about 3.4 A)."""
    a = np.array([r["C1'"].coord for r in structure[0]["A"]])
    b = np.array([r["C1'"].coord for r in structure[0]["B"]])[::-1]   # antiparallel pairing
    mid = (a + b) / 2
    return float(np.linalg.norm(mid[-1] - mid[0]) / (len(mid) - 1))


def pareto_front(sub, x, y):
    """Rows of `sub` not dominated on the two columns x and y (higher is better)."""
    keep = []
    for i, r in sub.iterrows():
        dominated = ((sub[x] >= r[x]) & (sub[y] >= r[y]) &
                     ((sub[x] > r[x]) | (sub[y] > r[y]))).any()
        keep.append(not dominated)
    return pd.Series(keep, index=sub.index)


def main(out_dir=RESULTS):
    meta = pd.read_csv(SUMMARY).set_index("probe")
    u = up_directions(N_ORIENTATIONS, SEED)

    rows = []
    for pdb in sorted(COMPLEX_DIR.glob("graphene_*.pdb")):
        probe = pdb.stem.split("_")[1]
        h, is_p = complex_heights(pdb)
        row = {"probe": probe,
               "gene": meta.at[probe, "gene"],
               "atoms": len(h),
               "phosphates": int(is_p.sum()),
               "height_above_sheet": round(float(h.max()), 2),
               "confidence": meta.at[probe, "confidence"],
               "plddt": meta.at[probe, "plddt"],
               "boltz_class": meta.at[probe, "boltz_class"]}

        # 1. rigid placement
        for name, lam in BUFFERS.items():
            row[f"f_{name}"] = round(100.0 * float((h <= lam).mean()), 1)
            row[f"P_{name}"] = int((h[is_p] <= lam).sum())

        # 2. probe anchored by its 3' end
        xyz, p_mask, _, terminal, _ = load_model(next(PROBE_DIR.glob(f"{probe}_*.cif")))
        H = anchored_heights(xyz, terminal, u)
        row["orientations"] = H.shape[1]
        for name, lam in BUFFERS.items():
            inside = H <= lam
            n_p = inside[p_mask].sum(axis=0)
            row[f"f_anchored_{name}"] = round(100.0 * float(np.median(inside.mean(axis=0))), 1)
            row[f"P_anchored_{name}"] = float(np.median(n_p))
            row[f"P_anchored_max_{name}"] = int(n_p.max())

        # 3. duplex anchored by the probe's 3' end: target phosphates inside lambda_D
        dup = next(DUPLEX_DIR.glob(f"{probe}_*.cif"), None) if DUPLEX_DIR.exists() else None
        if dup is not None:
            xyz, p_mask, chain, terminal, s = load_model(dup)
            H = anchored_heights(xyz, terminal, u)
            target_p = p_mask & (chain == "B")
            row["duplex_orientations"] = H.shape[1]
            row["duplex_rise_A"] = round(rise_per_bp(s), 2)
            for name, lam in BUFFERS.items():
                n_t = (H[target_p] <= lam).sum(axis=0)
                row[f"P_target_{name}"] = float(np.median(n_t))
                row[f"P_target_max_{name}"] = int(n_t.max())
        rows.append(row)

    cov = pd.DataFrame(rows).sort_values(["gene", "probe"])

    # final selection: Pareto front per gene of 3D confidence x anchored charge inside lambda_D
    key = "P_anchored_1x PBS"
    cov["final_shortlist"] = pd.concat([pareto_front(sub, "confidence", key)
                                        for _, sub in cov.groupby("gene")])

    # ";" so that the tables open in columns in Excel with a Portuguese locale
    cov.to_csv(out_dir / "debye_coverage.csv", index=False, sep=";")

    # aggregate by Boltz-2 confidence class
    agg = []
    for cls in ["GOOD", "MODERATE", "LOW"]:
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
            rec[f"{name} P median"] = sub[f"P_{name}"].median()
            rec[f"{name} P anchored median"] = sub[f"P_anchored_{name}"].median()
            if f"P_target_{name}" in sub:
                rec[f"{name} P target median"] = sub[f"P_target_{name}"].median()
        agg.append(rec)
    agg = pd.DataFrame(agg)
    agg.to_csv(out_dir / "debye_by_class.csv", index=False, sep=";")

    show = ["probe", "gene", "boltz_class", "confidence", "orientations",
            "f_1x PBS", "P_1x PBS", "P_anchored_1x PBS", "P_anchored_max_1x PBS"]
    show += [c for c in ("duplex_rise_A", "P_target_1x PBS", "P_target_0.1x PBS") if c in cov]
    print(cov[show + ["final_shortlist"]].to_string(index=False))
    print()
    print(agg.to_string(index=False))
    print()
    for name in BUFFERS:
        r1 = cov[f"f_{name}"].corr(cov["confidence"])
        r2 = cov[f"P_anchored_{name}"].corr(cov["confidence"])
        print(f"  Pearson r with Boltz confidence at {name:9s}: rigid atoms {r1:+.2f}, "
              f"anchored phosphates {r2:+.2f}")
    print(f"\n  final shortlist (Pareto, confidence x {key}): "
          f"{int(cov['final_shortlist'].sum())} of {len(cov)} probes")


if __name__ == "__main__":
    main()
