"""

Fase 0: Aquisição de Sequências, Alinhamento e Seleção de Regiões Conservadas
para alimentar o ViruScope e o pipeline GFET de design de probes.

"""

import os, sys, time, subprocess, shutil, platform
from pathlib import Path
from io import StringIO
from collections import Counter
from Bio import Entrez, SeqIO

# ── CONFIGURAÇÃO ─────────────────────────────────────────────────────────────
Entrez.email = "pg45861@uminho.pt"
Entrez.tool  = "gfet_probe_pipeline_v2"

OUTPUT_DIR   = Path("alignments")
MAX_SEQS     = 20
MIN_SEQ_LEN  = 200
MAX_SEQ_LEN  = 5000

# ── TARGETS ──────────────────────────────────────────────────────────────────
TARGETS = {
    "nuc": {
        "organism": "Staphylococcus aureus",
        "query_main": 'nuc[Gene Name] AND "Staphylococcus aureus"[Organism] AND 200:3000[Sequence Length]',
        "query_alt":  '"Staphylococcus aureus"[Organism] AND nuc[Title] AND 200:2000[Sequence Length]',
        "ref_accession": "NC_002952", "group": "A",
        "cross_react": ["S. epidermidis", "S. haemolyticus"],
        "notes": "Termostável nuclease. Benchmark GFET (Purwidyantri 2021).",
    },
            "rmpM": {
        "organism": "Neisseria meningitidis",
        "query_main": 'rmpM[Gene Name] AND "Neisseria meningitidis"[Organism] AND 100:1200[Sequence Length]',
        "query_alt": '("Neisseria meningitidis"[Organism]) AND (rmpM[Title] OR "class 4 outer membrane protein"[Title]) AND 100:1200[Sequence Length]',
        "ref_accession": "NC_003112", "group": "A",
        "cross_react": ["N. gonorrhoeae", "N. lactamica"],
        "notes": "OMP classe 4. Benchmark eletroquímico (Appaturi 2013).",
        "min_len": 100,
        "max_len": 1200,
    },
    "lytA": {
    "organism": "Streptococcus pneumoniae",
    "query_main": 'lytA[Gene Name] AND "Streptococcus pneumoniae"[Organism] AND 700:1300[Sequence Length]',
    "query_alt": '"Streptococcus pneumoniae"[Organism] AND lytA[Title] AND 700:1300[Sequence Length]',
    "ref_accession": "NC_003098", "group": "B",
    "cross_react": ["S. mitis", "S. oralis", "S. pseudopneumoniae"],
    "notes": "Autolisina LytA. Gold-standard S. pneumoniae.",
    "min_len": 700,
    "max_len": 1300,
    "cons_threshold": 0.80,
    "gap_threshold": 0.25,
    },
    "oprL": {
        "organism": "Pseudomonas aeruginosa",
        "query_main": 'oprL[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 300:2000[Sequence Length]',
        "query_alt":  '"Pseudomonas aeruginosa"[Organism] AND oprL[Title] AND 300:2000[Sequence Length]',
        "ref_accession": "NC_002516", "group": "B",
        "cross_react": ["P. fluorescens", "P. putida", "P. stutzeri"],
        "notes": "Lipoproteína OprL. Gold-standard P. aeruginosa.",
    },
    "algD": {
        "organism": "Pseudomonas aeruginosa",
        "query_main": 'algD[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 500:2500[Sequence Length]',
        "query_alt":  '"Pseudomonas aeruginosa"[Organism] AND algD[Title] AND 500:2500[Sequence Length]',
        "ref_accession": "NC_002516", "group": "B",
        "cross_react": ["P. fluorescens (algD homólogo)"],
        "notes": "AlgD. Marcador de fenótipo mucoide/biofilm.",
    },
    "frdB": {
        "organism": "Haemophilus influenzae",
        "query_main": 'frdB[Gene Name] AND "Haemophilus influenzae"[Organism] AND 200:2000[Sequence Length]',
        "query_alt":  '"Haemophilus influenzae"[Organism] AND frdB[Title] AND 200:2000[Sequence Length]',
        "ref_accession": "NC_000907", "group": "B",
        "cross_react": ["H. parainfluenzae", "H. haemolyticus"],
        "notes": "Fumarato redutase FrdB. Maior gap da literatura.",
    },
}

# ── ENCONTRAR MAFFT ─────────────────────────────────────────
MAFFT_CMD = r"MAFFT\mafft-7.526-win64-signed\mafft-win\mafft.bat"

def check_mafft():
    if Path(MAFFT_CMD).is_file():
        return True
    print("MAFFT não encontrado.")
    return False

# ── FETCH NCBI ────────────────────────────────────────────────────────────────
def safe_fetch(query, retmax=MAX_SEQS):
    """Descarrega sequências do NCBI com retry."""
    for attempt in range(3):
        try:
            time.sleep(0.4)
            handle  = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
            record  = Entrez.read(handle); handle.close()
            ids = record.get("IdList", [])
            if not ids:
                return None
            time.sleep(0.4)
            handle = Entrez.efetch(db="nucleotide", id=",".join(ids),
                                   rettype="fasta", retmode="text")
            data = handle.read(); handle.close()
            return data
        except Exception as e:
            print(f"    ⚠ Tentativa {attempt+1}/3 falhou: {e}")
            time.sleep(2 * (attempt + 1))
    return None

def filter_sequences(fasta_text, gene_name, min_len=MIN_SEQ_LEN, max_len=MAX_SEQ_LEN):
    lines = fasta_text.splitlines()

    start = 0
    while start < len(lines) and not lines[start].startswith((">", "#", "!", ";")):
        start += 1

    cleaned = "\n".join(lines[start:]).strip()
    if not cleaned:
        return []

    records, seen = [], set()
    for rec in SeqIO.parse(StringIO(cleaned), "fasta-blast"):
        if rec.id in seen:
            continue
        if not (min_len <= len(rec.seq) <= max_len):
            continue
        seen.add(rec.id)
        rec.id = f"{rec.id.split('.')[0]}_{gene_name}"
        rec.description = ""
        records.append(rec)
    return records

# ── ALINHAMENTO ───────────────────────────────────────────────────────────────
def run_mafft(input_fasta: Path, output_fasta: Path, n_seqs: int):
    """Corre MAFFT local. Devolve True se bem sucedido."""
    algo = "--auto" if n_seqs <= 200 else "--retree 1"
    cmd  = [MAFFT_CMD, algo, "--thread", "4",
            "--preservecase", str(input_fasta)]
    print(f"    ▶ {MAFFT_CMD} {algo} {input_fasta.name} → {output_fasta.name}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"    ✗ MAFFT retornou erro:\n{result.stderr[:400]}")
            return False
        output_fasta.write_text(result.stdout, encoding="utf-8")
        return True
    except FileNotFoundError:
        print(f"    ✗ Executável '{MAFFT_CMD}' não encontrado.")
        check_mafft()
        return False
    except subprocess.TimeoutExpired:
        print(f"    ✗ MAFFT excedeu o tempo limite (300s). Tentar com menos sequências.")
        return False

# ── CONSERVAÇÃO + JANELAS ─────────────────────────────────────────────────────
def analyse_conservation(aligned_fasta: Path, gene_name: str, target_info: dict):
    """Calcula conservação por posição e identifica janelas candidatas 20-25 nt."""
    records = list(SeqIO.parse(aligned_fasta, "fasta"))
    if len(records) < 2:
        return [], []
    n_seqs  = len(records)
    aln_len = len(records[0].seq)

    conservation = []
    for i in range(aln_len):
        col       = [str(r.seq[i]).upper() for r in records]
        top_count = Counter(col).most_common(1)[0][1]
        gap_freq  = col.count("-") / n_seqs
        conservation.append((top_count / n_seqs, gap_freq))

    # Sliding window
    all_windows = []
    for start in range(aln_len):
        for length in range(20, 26):
            end = start + length
            if end > aln_len:
                break
            window   = conservation[start:end]
            avg_cons = sum(c for c, g in window) / length
            avg_gap  = sum(g for c, g in window) / length
            cons_thr = target_info.get("cons_threshold", 0.85)
            gap_thr = target_info.get("gap_threshold", 0.20)

            if avg_cons >= cons_thr and avg_gap <= gap_thr:
                consensus = ""
                for i in range(start, end):
                    col = [str(r.seq[i]).upper() for r in records
                           if str(r.seq[i]) != "-"]
                    consensus += Counter(col).most_common(1)[0][0] if col else "N"
                gc = (consensus.count("G") + consensus.count("C")) / len(consensus)
                all_windows.append({
                    "start": start, "end": end, "length": length,
                    "avg_conservation": round(avg_cons, 3),
                    "avg_gap": round(avg_gap, 3),
                    "gc_content": round(gc, 3),
                    "consensus": consensus,
                })

    # Ordenar e remover sobreposições
    all_windows.sort(key=lambda x: -x["avg_conservation"])
    best, used = [], set()
    for w in all_windows:
        overlap = len(set(range(w["start"], w["end"])) & used)
        if overlap < w["length"] // 2:
            best.append(w)
            used.update(range(w["start"], w["end"]))
        if len(best) >= 10:
            break
    return conservation, best

# ── VIROSCOPE INPUT ───────────────────────────────────────────────────────────
def write_viroscope_fasta(windows: list, gene_name: str, gene_dir: Path):
    """
    Escreve um ficheiro FASTA com as sequências consenso das janelas candidatas,
    pronto para usar como input de probes no ViruScope.
    Formato esperado pelo ViruScope: >ID descrição \\n SEQUENCE
    """
    vs_path = gene_dir / f"{gene_name}_viroscope_probes.fasta"
    with open(vs_path, "w") as f:
        for i, w in enumerate(windows, 1):
            f.write(f">{gene_name}_probe_{i:02d}_pos{w['start']}-{w['end']}_"
                    f"cons{w['avg_conservation']}_gc{w['gc_content']}\n")
            f.write(f"{w['consensus']}\n")
    return vs_path

# ── RELATÓRIOS ────────────────────────────────────────────────────────────────
def write_outputs(gene_dir, gene_name, target_info, n_seqs, conservation, windows):
    # ── conserved_report.txt ─────────────────────────────────────────────
    rpt = gene_dir / "conserved_report.txt"
    with open(rpt, "w", encoding="utf-8") as f:
        f.write(f"{'═'*64}\n")
        f.write(f"  RELATÓRIO: {gene_name.upper()}\n")
        f.write(f"  {target_info['organism']}  |  Grupo {target_info['group']}\n")
        f.write(f"  {target_info['notes']}\n")
        f.write(f"{'═'*64}\n\n")
        f.write(f"Sequências alinhadas : {n_seqs}\n")
        f.write(f"Comprimento alinhamento: {len(conservation)} posições\n")
        if conservation:
            avg_g = sum(c for c,g in conservation) / len(conservation)
            h90   = sum(1 for c,g in conservation if c >= 0.90)
            f.write(f"Conservação média global : {avg_g:.3f}\n")
            f.write(f"Posições ≥ 90% conservadas: {h90}/{len(conservation)} ({100*h90/len(conservation):.1f}%)\n\n")
        f.write(f"{'─'*64}\n")
        f.write(f"  TOP {len(windows)} JANELAS CANDIDATAS PARA PROBE (20-25 nt)\n")
        cons_thr = target_info.get("cons_threshold", 0.85)
        gap_thr = target_info.get("gap_threshold", 0.20)
        print(f" ✔ {len(windows)} janelas candidatas (cons ≥ {cons_thr:.2f}, gaps ≤ {gap_thr:.0%})")

        f.write(f" Critérios: conservação ≥ {cons_thr:.2f}, gaps ≤ {gap_thr:.0%}, GC 40-60%\n")
        f.write(f"{'─'*64}\n\n")
        for i, w in enumerate(windows, 1):
            gc_flag = "⚠ GC fora de 40-60%" if not (0.40 <= w["gc_content"] <= 0.60) else "✔ GC ok"
            f.write(f"  #{i:02d}  pos {w['start']:>5}–{w['end']:>5}  len={w['length']} nt  "
                    f"cons={w['avg_conservation']:.3f}  gap={w['avg_gap']:.3f}  "
                    f"GC={w['gc_content']:.0%}  {gc_flag}\n")
            f.write(f"       5'-{w['consensus']}-3'\n\n")
        f.write(f"{'─'*64}\n")
        f.write(f"  REATIVIDADE CRUZADA — BLAST obrigatório contra:\n")
        for sp in target_info.get("cross_react", []):
            f.write(f"    • {sp}\n")
        f.write(f"\n  URL: https://blast.ncbi.nlm.nih.gov/Blast.cgi\n")
        f.write(f"  Query: cada sequência consenso acima\n")
        f.write(f"  Database: nr | Organism: excluir {target_info['organism']}\n\n")
        f.write(f"{'─'*64}\n")
        f.write(f"  PRÓXIMOS PASSOS\n")
        f.write(f"{'─'*64}\n")
        f.write(f"  1. Abrir aligned.fasta no Geneious/AliView para validação visual\n")
        f.write(f"  2. Importar {gene_name}_viroscope_probes.fasta no ViruScope\n")
        f.write(f"     → ViruScope faz scoring, anotação e validação in-silico\n")
        f.write(f"  3. Para cada janela seleccionada: NUPACK online (nupack.org)\n")
        f.write(f"     → ΔG < -12 kcal/mol, ensemble defect < 0.10\n")
        f.write(f"  4. BLAST de especificidade (URL acima)\n")
        f.write(f"  5. Probe final → Fase 1 workflow (3dDNA → AMBER → docking)\n")

    # ── probe_windows.tsv ────────────────────────────────────────────────
    tsv = gene_dir / "probe_windows.tsv"
    with open(tsv, "w", encoding="utf-8") as f:
        f.write("rank\tgene\tstart\tend\tlength\tavg_conservation\t"
                "avg_gap\tgc_content\tconsensus_5to3\n")
        for i, w in enumerate(windows, 1):
            f.write(f"{i}\t{gene_name}\t{w['start']}\t{w['end']}\t{w['length']}\t"
                    f"{w['avg_conservation']}\t{w['avg_gap']}\t{w['gc_content']}\t"
                    f"{w['consensus']}\n")

    # ── viroscope input ───────────────────────────────────────────────────
    vs_path = write_viroscope_fasta(windows, gene_name, gene_dir)

    return rpt, tsv, vs_path

# ── PIPELINE POR TARGET ────────────────────────────────────────────────────────
def target(gene_name, target_info, skip_fetch=False):
    print(f"\n{'━'*60}")
    print(f"  TARGET: {gene_name}  ({target_info['organism']}  |  Grupo {target_info['group']})")
    print(f"{'━'*60}")

    gene_dir = OUTPUT_DIR / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    raw_fasta     = gene_dir / "raw_sequences.fasta"
    aligned_fasta = gene_dir / "aligned.fasta"

    # 1. Descarregar (ou usar FASTA já existente)
    print(f"  [1/4] Descarregando sequências do NCBI...")
    if skip_fetch:
        if not raw_fasta.exists():
            print(f"    ✗ Ficheiro não encontrado: {raw_fasta}")
            print(f"    → Descarregar manualmente do NCBI e guardar em {raw_fasta}")
            return False
        print(f"    ℹ --skip-fetch: a usar {raw_fasta} já existente")
        records = list(SeqIO.parse(raw_fasta, "fasta"))
        print(f"    ✔ {len(records)} sequências carregadas")
    else:
        data = safe_fetch(target_info["query_main"])
        if not data:
            print(f"    ⚠ Query principal sem resultados. A tentar query alternativa...")
            data = safe_fetch(target_info["query_alt"])
        if not data:
            print(f"    ✗ NCBI não devolveu resultados para {gene_name}.")
            print(f"    → Descarregar manualmente em https://www.ncbi.nlm.nih.gov/nucleotide")
            print(f"    → Query sugerida: {target_info['query_alt']}")
            print(f"    → Guardar como {raw_fasta} e re-correr com --skip-fetch")
            return False

        records = filter_sequences(
            data,
            gene_name,
            min_len=target_info.get("min_len", MIN_SEQ_LEN),
            max_len=target_info.get("max_len", MAX_SEQ_LEN),
        )
        print(f"    ✔ {len(records)} sequências após filtragem ({target_info.get('min_len', MIN_SEQ_LEN)}–{target_info.get('max_len', MAX_SEQ_LEN)} bp)")

        if len(records) < 2:
            print("    ✗ Menos de 2 sequências. Alinhamento impossível.")
            return False

        SeqIO.write(records, raw_fasta, "fasta")
    # 2. Alinhar
    print(f"  [2/4] Alinhamento MAFFT ({len(records)} sequências)...")
    if not run_mafft(raw_fasta, aligned_fasta, len(records)):
        return False
    aligned_records = list(SeqIO.parse(aligned_fasta, "fasta"))
    if not aligned_records:
        print("    ✗ O MAFFT não produziu sequências alinhadas válidas.")
        return False
    print(f"    ✔ Alinhamento: {len(aligned_records)} seqs × {len(aligned_records[0].seq)} posições")

    # 3. Conservação
    print(f"  [3/4] Análise de conservação...")
    conservation, windows = analyse_conservation(aligned_fasta, gene_name, target_info)
    cons_thr = target_info.get("cons_threshold", 0.85)
    gap_thr = target_info.get("gap_threshold", 0.20)
    print(f"    ✔ {len(windows)} janelas candidatas (cons ≥ {cons_thr:.2f}, gaps ≤ {gap_thr:.0%})")

    # 4. Outputs
    print(f"  [4/4] A escrever outputs...")
    rpt, tsv, vs = write_outputs(
        gene_dir, gene_name, target_info, len(records), conservation, windows)
    print(f"    ✔ {rpt.name}")
    print(f"    ✔ {tsv.name}")
    print(f"    ✔ {vs.name}  ← importar no ViruScope")
    print(f"    ✔ {aligned_fasta.name}  ← importar no Geneious/AliView")
    return True

# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    # ── flags especiais ────────────────────────────────────────────────────
    skip_fetch = "--skip-fetch" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]

    # Verificar MAFFT antes de começar
    if not check_mafft():
        sys.exit(1)
    print(f"  ✔ MAFFT encontrado: {MAFFT_CMD}")
    if skip_fetch:
        print(f"  ℹ Modo --skip-fetch: a usar FASTA já existentes em alignments/<gene>/raw_sequences.fasta")

    targets_to_run = args if args else list(TARGETS.keys())
    invalid = [t for t in targets_to_run if t not in TARGETS]
    if invalid:
        print(f"✗ Targets desconhecidos: {invalid}")
        print(f"  Disponíveis: {list(TARGETS.keys())}")
        sys.exit(1)

    OUTPUT_DIR.mkdir(exist_ok=True)
    print(f"\n{'═'*60}")
    print(f"  GFET Probe Pipeline v2 — Fase 0")
    print(f"  OS: {platform.system()} | MAFFT: {MAFFT_CMD}")
    print(f"  Targets: {targets_to_run}")
    print(f"{'═'*60}")

    results = {}
    for gene in targets_to_run:
        results[gene] = target(gene, TARGETS[gene], skip_fetch=skip_fetch)

    print(f"\n{'═'*60}")
    print("  RESUMO FINAL")
    print(f"{'═'*60}")
    for gene, ok in results.items():
        status = "✔ OK    " if ok else "✗ FALHOU"
        print(f"  {status}  {gene:10s}  ({TARGETS[gene]['organism']})")

    vs_files = list(OUTPUT_DIR.glob("**/*_viroscope_probes.fasta"))
    if vs_files:
        print(f"\n  Ficheiros prontos para ViruScope:")
        for f in vs_files:
            print(f"    → {f}")
    print(f"\n  Todos os outputs em: {OUTPUT_DIR.resolve()}/")

if __name__ == "__main__":
    main()
