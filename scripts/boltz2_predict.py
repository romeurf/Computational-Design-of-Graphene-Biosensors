"""
boltz2_predict.py
─────────────────────────────────────────────────────────────────────────────
Script standalone: corre o Boltz-2 localmente para gerar estruturas 3D de
probes DNA e calcula o RMSD entre réplicas como critério de qualidade.

Uso:
    python boltz2_predict.py [--input CSV] [--filter {basic,nupack,all}]
                             [--samples N] [--device {cpu,cuda}]
                             [--out-csv boltz2_results.csv]

Argumentos:
    --input   : CSV com probes (padrão: alignments/FINAL_PROBES_ALL.csv)
    --filter  : quais probes processar
                  basic  → pass_basic == True   (padrão)
                  nupack → pass_nupack == True
                  all    → todas as linhas com sequência
    --samples : número de réplicas Boltz-2 por probe (padrão: 3)
    --device  : cpu ou cuda (padrão: cpu)
    --out-csv : CSV de resultados (padrão: boltz2_results.csv)
    --rmsd-max: RMSD máximo para pass_3d (padrão: 2.0 Å)

Saídas:
    structures/<probe_id>/boltz/*.cif     — estruturas geradas
    <out-csv>                             — CSV com probe_id, estruturas,
                                            rmsd_intra, pass_3d, status

Critério de qualidade (supervisor):
    Só modelos com variação pequena (RMSD < rmsd_max) avançam para
    Fase 1: 3dDNA → AMBER → docking molecular.

Referências:
    Wohlwend et al. 2024 — Boltz-2 (github.com/jwohlwend/boltz)
    Hamelryck & Manderick 2003 — Bio.PDB Superimposer
─────────────────────────────────────────────────────────────────────────────
"""

import argparse, csv, itertools, shutil, subprocess, sys
from pathlib import Path
from typing import Optional

import yaml

BASE_DIR   = Path(__file__).parent.parent          # raiz do repo
STRUCT_DIR = BASE_DIR / "output" / "structures"
STRUCT_DIR.mkdir(parents=True, exist_ok=True)

_RMSD_MAX_DEFAULT = 2.0  # Å


def calc_rmsd(cif_paths: list[Path]) -> Optional[float]:
    """
    RMSD médio entre todos os pares de estruturas CIF usando Bio.PDB Superimposer.
    Usa apenas os átomos C4' (esqueleto de desoxirribose).
    Devolve None se < 2 estruturas ou < 3 átomos em comum.
    Ref: Hamelryck & Manderick 2003 (Bio.PDB)
    """
    if len(cif_paths) < 2:
        return None
    try:
        from Bio.PDB import PDBParser, Superimposer
    except ImportError:
        print("    ⚠ Biopython não disponível — RMSD não calculado.")
        return None

    parser = PDBParser(QUIET=True)
    sup    = Superimposer()
    rmsds  = []
    for a, b in itertools.combinations(cif_paths, 2):
        try:
            s1 = parser.get_structure("s1", str(a))
            s2 = parser.get_structure("s2", str(b))
            at1 = [x for x in s1.get_atoms() if x.get_name() == "C4'"]
            at2 = [x for x in s2.get_atoms() if x.get_name() == "C4'"]
            n = min(len(at1), len(at2))
            if n < 3:
                continue
            sup.set_atoms(at1[:n], at2[:n])
            rmsds.append(sup.rms)
        except Exception:
            continue
    return round(sum(rmsds) / len(rmsds), 3) if rmsds else None


def _boltz_cmd() -> list[str] | None:
    if shutil.which("boltz"):
        return ["boltz"]
    try:
        import importlib.util
        if importlib.util.find_spec("boltz") is not None:
            return [sys.executable, "-m", "boltz"]
    except Exception:
        pass
    return None

def run_boltz(seq: str, probe_id: str, samples: int, device: str) -> list[Path]:
    """
    Cria o YAML de input e invoca `boltz predict`.
    Devolve lista de CIFs gerados (pode ser vazia se boltz falhar).
    """
    cmd_base = _boltz_cmd()
    if cmd_base is None:
        print("    boltz nao encontrado. Instalar: pip install boltz")
        return []

    out_dir = STRUCT_DIR / probe_id / "boltz"
    out_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = out_dir / "input.yaml"
    yaml_path.write_text(yaml.dump({
        "version": 1,
        "sequences": [{"dna": {"id": ["A"], "sequence": seq}}],
    }))

    cmd = cmd_base + [
        "predict", str(yaml_path),
        "--out_dir",          str(out_dir),
        "--recycling_steps",  "3",
        "--sampling_steps",   "200",
        "--diffusion_samples", str(samples),
        "--accelerator",      device,
        "--model",            "boltz2",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=900)
    except subprocess.CalledProcessError as e:
        print(f"    ✘ boltz retornou erro (probe {probe_id}):\n{e.stderr[-300:]}")
        return []
    except subprocess.TimeoutExpired:
        print(f"    ✘ boltz timeout (probe {probe_id}).")
        return []
    except Exception as e:
        print(f"    ✘ boltz erro inesperado: {e}")
        return []

    cifs = sorted(out_dir.glob("**/*.cif"))
    return cifs


def load_probes(csv_path: Path, filter_mode: str) -> list[dict]:
    if not csv_path.exists():
        sys.exit(f"✘ Ficheiro não encontrado: {csv_path}")
    probes = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            seq = row.get("sequence", "").strip().upper()
            if not seq:
                continue
            if filter_mode == "basic"  and str(row.get("pass_basic",  "")).lower() != "true":
                continue
            if filter_mode == "nupack" and str(row.get("pass_nupack", "")).lower() != "true":
                continue
            probes.append(row)
    return probes


def run(input_csv: Path, filter_mode: str, samples: int,
        device: str, out_csv: Path, rmsd_max: float):

    probes = load_probes(input_csv, filter_mode)
    print(f"\n{'═'*60}")
    print(f"  Boltz-2 Predict — {len(probes)} probes | {samples} réplicas | device={device}")
    print(f"  RMSD threshold: {rmsd_max} Å   filtro: {filter_mode}")
    print(f"{'═'*60}")

    results = []
    for i, probe in enumerate(probes, 1):
        pid = probe.get("probe_id", f"probe{i:03d}")
        seq = probe["sequence"]
        print(f"\n  [{i}/{len(probes)}] {pid}  ({len(seq)} nt)")

        # Skip se já temos CIFs suficientes
        existing_cifs = sorted((STRUCT_DIR / pid / "boltz").glob("**/*.cif")) \
                        if (STRUCT_DIR / pid / "boltz").exists() else []
        if len(existing_cifs) >= samples:
            print(f"    ℹ {len(existing_cifs)} CIFs já existem — a re-calcular RMSD.")
            cifs = existing_cifs
        else:
            cifs = run_boltz(seq, pid, samples, device)

        if not cifs:
            results.append({**probe,
                             "structure_cif": "",
                             "boltz_cifs":    "",
                             "n_models":      0,
                             "rmsd_intra":    "",
                             "pass_3d":       False,
                             "status":        "failed",
                             "notes":         "boltz falhou ou não instalado"})
            continue

        print(f"    ✔ {len(cifs)} modelos gerados")

        rmsd = calc_rmsd(cifs)
        pass_3d = (rmsd is not None and rmsd < rmsd_max)
        rmsd_str = f"{rmsd:.3f} Å" if rmsd is not None else "N/A"

        if rmsd is None:
            flag = "⚠ RMSD não calculado"
        elif pass_3d:
            flag = f"✔ PASS (RMSD={rmsd_str})"
        else:
            flag = f"✘ FAIL (RMSD={rmsd_str} ≥ {rmsd_max} Å)"
        print(f"    {flag}")

        results.append({**probe,
                         "structure_cif": str(cifs[0]),
                         "boltz_cifs":    ";".join(str(c) for c in cifs),
                         "n_models":      len(cifs),
                         "rmsd_intra":    rmsd if rmsd is not None else "",
                         "pass_3d":       pass_3d,
                         "status":        "done",
                         "notes":         f"RMSD={rmsd_str}"})

    if results:
        fields = list(results[0].keys())
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader(); w.writerows(results)

    n_pass  = sum(1 for r in results if r.get("pass_3d") is True)
    n_total = len(results)
    print(f"\n  ✔ Concluído: {n_pass}/{n_total} probes pass_3d  →  {out_csv}")
    print(f"  Critério: RMSD entre réplicas < {rmsd_max} Å")
    print(f"  Próximo passo para PASS: 3dDNA → AMBER → docking molecular")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Boltz-2 standalone prediction + RMSD filter")
    ap.add_argument("--input",    default="output/FINAL_PROBES_ALL.csv")
    ap.add_argument("--filter",   default="basic", choices=["basic", "nupack", "all"])
    ap.add_argument("--samples",  type=int, default=3,
                    help="Número de réplicas Boltz-2 (padrão: 3)")
    ap.add_argument("--device",   default="cpu", choices=["cpu", "cuda"])
    ap.add_argument("--out-csv",  default="output/boltz2_results.csv")
    ap.add_argument("--rmsd-max", type=float, default=_RMSD_MAX_DEFAULT,
                    help=f"RMSD máximo para pass_3d (padrão: {_RMSD_MAX_DEFAULT} Å)")
    args = ap.parse_args()
    run(
        input_csv   = BASE_DIR / args.input,
        filter_mode = args.filter,
        samples     = args.samples,
        device      = args.device,
        out_csv     = BASE_DIR / args.out_csv,
        rmsd_max    = args.rmsd_max,
    )
