"""
nucleofold_predict.py
─────────────────────────────────────────────────────────────────────────────
Script standalone: submete probes à API pública do NucleoFold (3dRNA)
e guarda as estruturas CIF em structures/<probe_id>_nucleofold.cif.

Uso:
    python nucleofold_predict.py [--input CSV] [--filter {basic,nupack,all}]
                                 [--out-csv nucleofold_results.csv]

Argumentos:
    --input   : CSV com probes (padrão: alignments/FINAL_PROBES_ALL.csv)
    --filter  : quais probes processar
                  basic  → pass_basic == True   (padrão)
                  nupack → pass_nupack == True
                  all    → todas as linhas com sequência
    --out-csv : onde guardar o CSV de resultados

Saídas:
    structures/<probe_id>_nucleofold.cif   — estruturas em formato CIF
    <out-csv>                              — CSV com colunas probe_id,
                                             structure_cif, status, notes

Referência do servidor:
    3dRNA v2.0 — http://biophy.hust.edu.cn/new/3dRNA
    Zhao et al. 2012, NAR; Wang et al. 2023 (update)
─────────────────────────────────────────────────────────────────────────────
"""

import argparse, csv, sys, time
from pathlib import Path
from typing import Optional

import requests

BASE_DIR   = Path(__file__).parent
STRUCT_DIR = BASE_DIR / "structures"
STRUCT_DIR.mkdir(parents=True, exist_ok=True)

# NucleoFold / 3dRNA API endpoints
_API_SUBMIT = "http://biophy.hust.edu.cn/new/3dRNA/api/submit"
_API_RESULT = "http://biophy.hust.edu.cn/new/3dRNA/api/result"
_POLL_SEC   = 10    # segundos entre polls
_MAX_POLLS  = 72    # 72 × 10 s = 12 min máx por probe


def submit_nucleofold(sequence: str, probe_id: str) -> Optional[str]:
    """Submete uma sequência DNA ao NucleoFold. Devolve job_id ou None."""
    try:
        resp = requests.post(
            _API_SUBMIT,
            json={"sequence": sequence, "type": "DNA", "name": probe_id},
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        return data.get("job_id") or data.get("id")
    except Exception as e:
        print(f"    ✘ Submissão falhou ({probe_id}): {e}")
        return None


def poll_nucleofold(job_id: str, probe_id: str) -> Optional[str]:
    """
    Aguarda conclusão do job NucleoFold e devolve o URL do ficheiro CIF,
    ou None se falhar ou expirar.
    """
    for i in range(_MAX_POLLS):
        time.sleep(_POLL_SEC)
        try:
            r = requests.get(f"{_API_RESULT}/{job_id}", timeout=30)
            r.raise_for_status()
            d = r.json()
            status = d.get("status", "").lower()
            if status in ("done", "finished", "completed"):
                return d.get("cif_url") or d.get("structure_url") or d.get("result_url")
            if status in ("error", "failed"):
                print(f"    ✘ Job {job_id} falhou (probe {probe_id}).")
                return None
            if i % 6 == 0:
                print(f"    ⏳ Aguardando NucleoFold ({probe_id}) — {(i+1)*_POLL_SEC}s ...")
        except Exception as e:
            print(f"    ⚠ Erro ao consultar job {job_id}: {e}")
    print(f"    ✘ Timeout NucleoFold ({probe_id}).")
    return None


def download_cif(url: str, probe_id: str) -> Optional[Path]:
    """Descarrega o CIF e guarda em structures/<probe_id>_nucleofold.cif."""
    out = STRUCT_DIR / f"{probe_id}_nucleofold.cif"
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        out.write_text(r.text)
        return out
    except Exception as e:
        print(f"    ✘ Download CIF falhou ({probe_id}): {e}")
        return None


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


def run(input_csv: Path, filter_mode: str, out_csv: Path):
    probes = load_probes(input_csv, filter_mode)
    print(f"\n{'═'*60}")
    print(f"  NucleoFold Predict — {len(probes)} probes  |  filtro: {filter_mode}")
    print(f"{'═'*60}")

    results = []
    for i, probe in enumerate(probes, 1):
        pid = probe.get("probe_id", f"probe{i:03d}")
        seq = probe["sequence"]
        print(f"\n  [{i}/{len(probes)}] {pid}  ({len(seq)} nt)")

        # Skip se já temos CIF
        existing = STRUCT_DIR / f"{pid}_nucleofold.cif"
        if existing.exists():
            print(f"    ℹ CIF já existe — a saltar.")
            results.append({**probe, "structure_cif": str(existing),
                             "status": "existing", "notes": ""})
            continue

        job_id = submit_nucleofold(seq, pid)
        if not job_id:
            results.append({**probe, "structure_cif": "",
                             "status": "submit_failed", "notes": "submit failed"})
            continue

        print(f"    ✔ job_id={job_id}  a aguardar resultado...")
        cif_url = poll_nucleofold(job_id, pid)
        if not cif_url:
            results.append({**probe, "structure_cif": "",
                             "status": "failed", "notes": f"job_id={job_id}"})
            continue

        cif_path = download_cif(cif_url, pid)
        if cif_path:
            print(f"    ✔ {cif_path.name}")
            results.append({**probe, "structure_cif": str(cif_path),
                             "status": "done", "notes": ""})
        else:
            results.append({**probe, "structure_cif": "",
                             "status": "download_failed",
                             "notes": f"cif_url={cif_url}"})

    # Escrever CSV de resultados
    if results:
        fields = list(results[0].keys())
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader(); w.writerows(results)

    done  = sum(1 for r in results if r["status"] in ("done", "existing"))
    total = len(results)
    print(f"\n  ✔ Concluído: {done}/{total} estruturas  →  {out_csv}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="NucleoFold standalone prediction")
    ap.add_argument("--input",   default="alignments/FINAL_PROBES_ALL.csv",
                    help="CSV de probes (padrão: alignments/FINAL_PROBES_ALL.csv)")
    ap.add_argument("--filter",  default="basic",
                    choices=["basic", "nupack", "all"],
                    help="Quais probes processar (padrão: basic)")
    ap.add_argument("--out-csv", default="nucleofold_results.csv",
                    help="CSV de resultados (padrão: nucleofold_results.csv)")
    args = ap.parse_args()
    run(
        input_csv   = BASE_DIR / args.input,
        filter_mode = args.filter,
        out_csv     = BASE_DIR / args.out_csv,
    )
