"""
nucleofold_predict.py — wrapper para NucleoFold3D.py (pipeline 3D local)
─────────────────────────────────────────────────────────────────────────────
Prepara o CSV de input para NucleoFold3D.py, invoca o pipeline,
recolhe os PDB resultantes e calcula o RMSD entre as 5 réplicas
como critério de qualidade interna.

REQUISITO: NucleoFold3D.py precisa de ferramentas Linux:
  • nsp       ~/nsp_final/cirRNA_and_DNA/nsp  (coarse-grained assembly)
  • RNAfold   ViennaRNA                       (estrutura secundária)
  • tleap, sander, ambpdb  AmberTools ≥ 22    (minimização AMBER)
  → Em Windows: correr dentro do WSL.
  → Instalar AmberTools: conda install -c conda-forge ambertools
  → Instalar ViennaRNA:  conda install -c bioconda viennarna

Uso:
    # Passo 1 — preparar CSV (corre em Python nativo, Windows ou Linux)
    python nucleofold_predict.py --filter basic --prepare-only

    # Passo 2 — correr NucleoFold3D (requer WSL/Linux)
    #   A flag --run invoca automaticamente; sem ela imprime o comando manual.
    python nucleofold_predict.py --filter basic [--run]

    # Passo 3 — recolher resultados e calcular RMSD
    python nucleofold_predict.py --collect

Argumentos:
    --input        CSV de probes (padrão: alignments/FINAL_PROBES_ALL.csv)
    --filter       basic | nupack | all  (padrão: nupack)
    --nucleofold   path para NucleoFold3D.py
    --results-dir  directório de resultados do NucleoFold3D
                   (padrão: ~/scripts/results_gfet_probes)
    --rmsd-max     RMSD máximo para pass_3d entre as 5 réplicas (padrão: 2.0 Å)
    --prepare-only só cria o CSV, não corre NucleoFold3D
    --run          invoca NucleoFold3D.py directamente (requer WSL/Linux)
    --collect      recolhe resultados existentes e calcula RMSD
    --out-csv      CSV de resultados (padrão: nucleofold_results.csv)

Ref. NucleoFold3D:  RNAfold → NSP → AMBER OL15 (DNA) / OL3 (RNA)
                    Mathews 2004 DNA params; Onufriev GB igb=5; saltcon=0.15
─────────────────────────────────────────────────────────────────────────────
"""

import argparse, csv, itertools, subprocess, sys
from pathlib import Path
from typing import Optional

BASE_DIR   = Path(__file__).parent.parent          # raiz do repo
STRUCT_DIR = BASE_DIR / "output" / "structures"
STRUCT_DIR.mkdir(parents=True, exist_ok=True)

_RMSD_MAX_DEFAULT    = 2.0
_NUCLEOFOLD3D_DEFAULT = (
    Path.home() / "OneDrive" / "Ambiente de Trabalho" / "Tese" / "NucleoFold" / "NucleoFold3D.py"
)


# ─── RMSD via Biopython (PDB, usa C4' do esqueleto) ──────────────────────────

def calc_rmsd_pdbs(pdb_paths: list[Path]) -> Optional[float]:
    """
    RMSD médio entre pares de PDB usando Bio.PDB Superimposer (C4').
    Idêntico ao calc_rmsd de boltz2_predict.py mas lê PDB em vez de CIF.
    Ref: Hamelryck & Manderick 2003 (Bio.PDB)
    """
    if len(pdb_paths) < 2:
        return None
    try:
        from Bio.PDB import PDBParser, Superimposer
    except ImportError:
        print("    ⚠ Biopython não disponível — instalar: pip install biopython")
        return None
    parser = PDBParser(QUIET=True)
    sup    = Superimposer()
    rmsds  = []
    for a, b in itertools.combinations(pdb_paths, 2):
        try:
            s1 = parser.get_structure("s1", str(a))
            s2 = parser.get_structure("s2", str(b))
            at1 = [x for x in s1.get_atoms() if x.get_name() == "C4'"]
            at2 = [x for x in s2.get_atoms() if x.get_name() == "C4'"]
            n   = min(len(at1), len(at2))
            if n < 3:
                continue
            sup.set_atoms(at1[:n], at2[:n])
            rmsds.append(sup.rms)
        except Exception:
            continue
    return round(sum(rmsds) / len(rmsds), 3) if rmsds else None


# ─── Helpers ─────────────────────────────────────────────────────────────────

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


def prepare_input_csv(probes: list[dict], out_csv: Path):
    """Cria o CSV no formato que NucleoFold3D.py espera: job_name,sequence"""
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["job_name", "sequence"])
        for p in probes:
            pid = p.get("probe_id", "probe")
            seq = p["sequence"]
            w.writerow([pid, seq])
    print(f"  ✔ CSV de input criado: {out_csv}  ({len(probes)} probes)")


def collect_results(probes: list[dict], results_dir: Path,
                    rmsd_max: float, out_csv: Path):
    """
    Percorre resultados do NucleoFold3D.py, calcula RMSD entre as 5 réplicas
    e gera o CSV de resultados.
    """
    results = []
    for probe in probes:
        pid       = probe.get("probe_id", "probe")
        job_dir   = results_dir / pid

        # NucleoFold3D pode adicionar sufixo _2, _3 se já existia pasta
        # tentar encontrar o directório mais recente
        candidates = sorted(results_dir.glob(f"{pid}*"), key=lambda p: p.stat().st_mtime,
                            reverse=True) if results_dir.exists() else []
        if candidates:
            job_dir = candidates[0]

        if not job_dir.exists():
            print(f"  ⚠ {pid}: resultados não encontrados em {job_dir}")
            results.append({**probe, "structure_pdb": "", "nucleofold_pdbs": "",
                             "n_models": 0, "rmsd_intra": "", "pass_3d": False,
                             "status": "not_found", "notes": ""})
            continue

        # Recolher os 5 PDB de final_structures/
        fs_dir = job_dir / "final_structures"
        pdbs   = sorted(fs_dir.glob(f"{pid}[1-5].pdb")) if fs_dir.exists() else []
        if not pdbs:
            # Tentar com nome do job real (pode ter sufixo)
            actual_name = job_dir.name
            pdbs        = sorted(fs_dir.glob(f"{actual_name}[1-5].pdb")) if fs_dir.exists() else []

        # Estrutura de menor energia
        min_pdb = next(job_dir.glob(f"{pid}*_minenerg.pdb"), None)
        if min_pdb is None:
            min_pdb = next(job_dir.glob(f"*_minenerg.pdb"), None)

        if not pdbs and min_pdb is None:
            print(f"  ⚠ {pid}: nenhum PDB encontrado em {job_dir}")
            results.append({**probe, "structure_pdb": "", "nucleofold_pdbs": "",
                             "n_models": 0, "rmsd_intra": "", "pass_3d": False,
                             "status": "no_pdb", "notes": ""})
            continue

        # Copiar min energy PDB para structures/
        dest_pdb = None
        if min_pdb:
            dest_pdb = STRUCT_DIR / f"{pid}_nucleofold_minenerg.pdb"
            import shutil
            shutil.copy(min_pdb, dest_pdb)

        # RMSD entre as 5 réplicas (qualidade interna)
        rmsd = calc_rmsd_pdbs(pdbs) if len(pdbs) >= 2 else None
        pass_3d = (rmsd is not None and rmsd < rmsd_max)
        rmsd_str = f"{rmsd:.3f} Å" if rmsd is not None else "N/A (< 2 modelos)"

        flag = ("✔ PASS" if pass_3d else "✘ FAIL") if rmsd is not None else "⚠ RMSD N/A"
        print(f"  {flag}  {pid}  {len(pdbs)} modelos  RMSD intra={rmsd_str}")

        results.append({
            **probe,
            "structure_pdb":   str(dest_pdb) if dest_pdb else "",
            "nucleofold_pdbs": ";".join(str(p) for p in pdbs),
            "n_models":        len(pdbs),
            "rmsd_intra":      rmsd if rmsd is not None else "",
            "pass_3d":         pass_3d,
            "status":          "done",
            "notes":           f"RMSD={rmsd_str}",
        })

    if results:
        fields = list(results[0].keys())
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(results)

    n_pass = sum(1 for r in results if r.get("pass_3d") is True)
    print(f"\n  ✔ Recolhidos: {n_pass}/{len(results)} pass_3d  →  {out_csv}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def run(input_csv: Path, filter_mode: str, nucleofold_script: Path,
        results_dir: Path, rmsd_max: float, out_csv: Path,
        prepare_only: bool, run_now: bool, collect: bool):

    probes = load_probes(input_csv, filter_mode)
    print(f"\n{'═'*60}")
    print(f"  NucleoFold3D Wrapper — {len(probes)} probes  |  filtro: {filter_mode}")
    print(f"{'═'*60}")

    # CSV input para NucleoFold3D.py
    input_for_nf3d = BASE_DIR / "output" / "nucleofold_input.csv"

    if not collect:
        prepare_input_csv(probes, input_for_nf3d)

        cmd_str = f"python \"{nucleofold_script}\" --csv \"{input_for_nf3d}\""

        if prepare_only:
            print(f"\n  CSV de input pronto. Correr manualmente (em WSL/Linux):")
            print(f"    {cmd_str}")
            return

        if run_now:
            print(f"\n  A invocar NucleoFold3D.py ...")
            print(f"    {cmd_str}")
            try:
                subprocess.run(
                    ["python", str(nucleofold_script), "--csv", str(input_for_nf3d)],
                    check=True, timeout=86400   # 24 h
                )
            except FileNotFoundError:
                print("  ✘ python/NucleoFold3D.py não encontrado.")
                print("    Em Windows, correr dentro do WSL:")
                print(f"    wsl python {nucleofold_script} --csv {input_for_nf3d}")
                return
            except subprocess.CalledProcessError as e:
                print(f"  ✘ NucleoFold3D.py falhou: {e}")
                return
        else:
            print(f"\n  ─────────────────────────────────────────────────────")
            print(f"  Próximo passo: correr NucleoFold3D.py (requer WSL/Linux)")
            print(f"  ─────────────────────────────────────────────────────")
            print(f"  Abrir WSL e correr:")
            print(f"    python {nucleofold_script} --csv {input_for_nf3d}")
            print(f"  Depois recolher resultados com:")
            print(f"    python nucleofold_predict.py --collect")
            print(f"  ─────────────────────────────────────────────────────")
            return

    # Recolher resultados
    print(f"\n  A recolher resultados de: {results_dir}")
    collect_results(probes, results_dir, rmsd_max, out_csv)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="NucleoFold3D wrapper + RMSD intra-réplicas")
    ap.add_argument("--input",       default="output/FINAL_PROBES_ALL.csv")
    ap.add_argument("--filter",      default="nupack", choices=["basic", "nupack", "all"])
    ap.add_argument("--nucleofold",  default=str(_NUCLEOFOLD3D_DEFAULT),
                    help="Path para NucleoFold3D.py")
    ap.add_argument("--results-dir", default=str(Path.home() / "scripts" / "results_nucleofold_input"),
                    help="Directório de resultados do NucleoFold3D.py")
    ap.add_argument("--rmsd-max",    type=float, default=_RMSD_MAX_DEFAULT)
    ap.add_argument("--out-csv",     default="output/nucleofold_results.csv")
    ap.add_argument("--prepare-only", action="store_true",
                    help="Só criar CSV de input, não correr NucleoFold3D")
    ap.add_argument("--run",         action="store_true",
                    help="Invocar NucleoFold3D.py directamente (requer WSL/Linux)")
    ap.add_argument("--collect",     action="store_true",
                    help="Recolher resultados existentes e calcular RMSD")
    args = ap.parse_args()
    run(
        input_csv        = BASE_DIR / args.input,
        filter_mode      = args.filter,
        nucleofold_script = Path(args.nucleofold),
        results_dir      = Path(args.results_dir),
        rmsd_max         = args.rmsd_max,
        out_csv          = BASE_DIR / args.out_csv,
        prepare_only     = args.prepare_only,
        run_now          = args.run,
        collect          = args.collect,
    )
