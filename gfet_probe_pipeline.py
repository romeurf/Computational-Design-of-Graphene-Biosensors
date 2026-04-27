"""
gfet_probe_pipeline.py (ViruScope integrado)
=======================================================
Fase 0: Aquisição → Alinhamento → Conservação → Scoring de Probes
Integra a lógica de scoring do ViruScope (primer3-py):
  - Melting Temperature (NN method, 50 mM Na+, 37°C)
  - GC content
  - Hairpin ΔG (critério de no-fold para ssDNA)
  - Homodimer ΔG
  - Conservation score (do alinhamento MAFFT)

"""

from unittest import result

import pandas as pd
import os, sys, time, subprocess, shutil, platform
from pathlib import Path
from io import StringIO
from collections import Counter
from Bio import Entrez, SeqIO
import primer3
HAS_PRIMER3 = True


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURAÇÃO
# ══════════════════════════════════════════════════════════════════════════════
Entrez.email = "pg45861@uminho.pt"
Entrez.tool  = "gfet_probe_pipeline_v3"

OUTPUT_DIR = Path("alignments")
MAX_SEQS   = 20

# Critérios ViruScope / GFET (podem ser ajustados por target)
PROBE_LEN_MIN   = 20        # nt mínimo da probe
PROBE_LEN_MAX   = 25        # nt máximo da probe
TM_MIN          = 55.0      # °C — Tm mínima aceitável
TM_MAX          = 72.0      # °C — Tm máxima aceitável
GC_MIN          = 0.40      # GC mínimo
GC_MAX          = 0.60      # GC máximo
HAIRPIN_DG_MAX  = -2.0      # kcal/mol — acima disto = sem fold significativo (critério no-fold)
HOMODIMER_DG_MAX= -6.0      # kcal/mol — acima disto = baixo risco de dimerização

# Primer3 thermodynamic parameters (equivalentes ao ViruScope)
P3_MV_CONC  = 50.0    # [Na+] mM
P3_DV_CONC  = 0.0     # [Mg2+] mM
P3_DNTP     = 0.0     # [dNTP] mM
P3_DNA_CONC = 250.0   # [oligo] nM (250 nM típico para biossensores)
P3_TEMP_C   = 37.0    # temperatura de referência (condições fisiológicas/GFET)

# ══════════════════════════════════════════════════════════════════════════════
# TARGETS — parâmetros POR TARGET
# ══════════════════════════════════════════════════════════════════════════════
TARGETS = {
    "nuc": {
        "organism"       : "Staphylococcus aureus",
        "group"          : "A",
        "notes"          : "Nuclease termoestável. Benchmark GFET (Purwidyantri 2021).",
        "cross_react"    : ["S. epidermidis", "S. haemolyticus"],
        "ref_accession"  : "NC_002952",
        "min_len"        : 200,   # comprimento mínimo das sequências a aceitar (bp)
        "max_len"        : 3000,
        "cons_threshold" : 0.85,
        "gap_threshold"  : 0.20,
        "query_main"     : 'nuc[Gene Name] AND "Staphylococcus aureus"[Organism] AND 200:3000[Sequence Length]',
        "query_alt"      : '"Staphylococcus aureus"[Organism] AND nuc[Title] AND 200:2000[Sequence Length]',
    },
    "rmpM": {
        "organism"       : "Neisseria meningitidis",
        "group"          : "A",
        "notes"          : "OMP classe 4. Benchmark eletroquímico (Appaturi 2013).",
        "cross_react"    : ["N. gonorrhoeae", "N. lactamica"],
        "ref_accession"  : "NC_003112",
        "min_len"        : 100,
        "max_len"        : 1200,
        "cons_threshold" : 0.80,   # limiar mais permissivo — região problemática
        "gap_threshold"  : 0.30,
        # rmpM tem poucos homólogos no NCBI — usar ref + expanded search
        "query_main"     : 'rmpM[Gene Name] AND "Neisseria meningitidis"[Organism] AND 100:1200[Sequence Length]',
        "query_alt"      : '("Neisseria meningitidis"[Organism]) AND (rmpM[Title] OR "class 4 outer membrane protein"[Title]) AND 100:1200[Sequence Length]',
        "query_fallback" : '"Neisseria meningitidis"[Organism] AND "outer membrane protein"[Title] AND 100:1200[Sequence Length]',
        "allow_single_seq": True,  # ← fallback: se só 1 seq, gerar probes da referência
    },
    "lytA": {
        "organism"       : "Streptococcus pneumoniae",
        "group"          : "B",
        "notes"          : "Autolisina LytA. Gold-standard S. pneumoniae.",
        "cross_react"    : ["S. mitis", "S. oralis", "S. pseudopneumoniae"],
        "ref_accession"  : "NC_003098",
        "min_len"        : 700,
        "max_len"        : 1300,
        "cons_threshold" : 0.70,
        "gap_threshold"  : 0.40,
        "query_main"     : 'lytA[Gene Name] AND "Streptococcus pneumoniae"[Organism] AND 700:1300[Sequence Length]',
        "query_alt"      : '"Streptococcus pneumoniae"[Organism] AND lytA[Title] AND 700:1300[Sequence Length]',
    },
    "oprL": {
        "organism"       : "Pseudomonas aeruginosa",
        "group"          : "B",
        "notes"          : "Lipoproteína OprL. Gold-standard P. aeruginosa.",
        "cross_react"    : ["P. fluorescens", "P. putida", "P. stutzeri"],
        "ref_accession"  : "NC_002516",
        "min_len"        : 300,
        "max_len"        : 2000,
        "cons_threshold" : 0.85,
        "gap_threshold"  : 0.20,
        "query_main"     : 'oprL[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 300:2000[Sequence Length]',
        "query_alt"      : '"Pseudomonas aeruginosa"[Organism] AND oprL[Title] AND 300:2000[Sequence Length]',
    },
    "algD": {
        "organism"       : "Pseudomonas aeruginosa",
        "group"          : "B",
        "notes"          : "AlgD. Marcador de fenótipo mucoide/biofilm.",
        "cross_react"    : ["P. fluorescens (algD homólogo)"],
        "ref_accession"  : "NC_002516",
        "min_len"        : 500,
        "max_len"        : 2500,
        "cons_threshold" : 0.85,
        "gap_threshold"  : 0.20,
        "query_main"     : 'algD[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 500:2500[Sequence Length]',
        "query_alt"      : '"Pseudomonas aeruginosa"[Organism] AND algD[Title] AND 500:2500[Sequence Length]',
    },
    "frdB": {
        "organism"       : "Haemophilus influenzae",
        "group"          : "B",
        "notes"          : "Fumarato redutase FrdB. Maior gap da literatura.",
        "cross_react"    : ["H. parainfluenzae", "H. haemolyticus"],
        "ref_accession"  : "NC_000907",
        "min_len"        : 200,
        "max_len"        : 2000,
        "cons_threshold" : 0.85,
        "gap_threshold"  : 0.20,
        "query_main"     : 'frdB[Gene Name] AND "Haemophilus influenzae"[Organism] AND 200:2000[Sequence Length]',
        "query_alt"      : '"Haemophilus influenzae"[Organism] AND frdB[Title] AND 200:2000[Sequence Length]',
    },
}

# ══════════════════════════════════════════════════════════════════════════════
# MAFFT
# ══════════════════════════════════════════════════════════════════════════════
MAFFT_CMD = r"MAFFT\mafft-7.526-win64-signed\mafft-win\mafft.bat"

# ══════════════════════════════════════════════════════════════════════════════
# VIROSCOPE SCORING (primer3-py)
# ══════════════════════════════════════════════════════════════════════════════
def calc_tm(seq: str) -> float:
    """Tm por Nearest-Neighbour (SantaLucia 1998), condições fisiológicas."""
    if not HAS_PRIMER3:
        # Fallback: Wallace rule básica
        return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))
    return primer3.calc_tm(
        seq,
        mv_conc=P3_MV_CONC,
        dv_conc=P3_DV_CONC,
        dntp_conc=P3_DNTP,
        dna_conc=P3_DNA_CONC,
    )

def calc_hairpin_dg(seq: str) -> float:
    try:
        result = primer3.calc_hairpin(
            seq,
            mv_conc=P3_MV_CONC,
            dv_conc=P3_DV_CONC,
            dntp_conc=P3_DNTP,
            dna_conc=P3_DNA_CONC,
            temp_c=P3_TEMP_C,
        )

        # se não há estrutura → NÃO inventar 0.0
        if not result.structure_found:
            return None

        dg = result.dg / 1000

        if dg > 50 or dg < -50:
            return None

        return round(dg, 2)

    except Exception:
        return None

def calc_homodimer_dg(seq: str) -> float:
    """Homodimer ΔG a 37°C."""
    if not HAS_PRIMER3:
        return None

    try:
        result = primer3.calc_homodimer(
            seq,
            mv_conc=P3_MV_CONC,
            dv_conc=P3_DV_CONC,
            dntp_conc=P3_DNTP,
            dna_conc=P3_DNA_CONC,
            temp_c=P3_TEMP_C,
        )

        # se não há estrutura → não inventar 0.0
        if not result.structure_found:
            return None

        dg = result.dg / 1000

        # sanity check (evita valores absurdos)
        if dg > 50 or dg < -50:
            return None

        return round(dg, 2)

    except Exception:
        return None

def score_probe(seq: str) -> dict:
    """
    Calcula todos os scores ViruScope-style para uma sequência de probe.
    Retorna dict com métricas e flags de PASS/FAIL.
    """
    seq = seq.upper().replace("U", "T")   # aceitar RNA também
    n   = len(seq)
    gc  = (seq.count("G") + seq.count("C")) / n
    tm  = calc_tm(seq)
    hp  = calc_hairpin_dg(seq)
    hd  = calc_homodimer_dg(seq)

    # Critérios de PASS (equivalentes ao ViruScope)
    gc_ok = GC_MIN <= gc <= GC_MAX
    tm_ok = TM_MIN <= tm <= TM_MAX
    fold_ok = (hp is None) or (hp > HAIRPIN_DG_MAX)       # ΔG > -2 → não faz fold significativo
    dimer_ok = hd > HOMODIMER_DG_MAX     # ΔG > -6 → baixo risco de dimerização

    pass_all = gc_ok and tm_ok and fold_ok and dimer_ok

    return {
        "length"        : n,
        "gc_content"    : round(gc, 3),
        "tm_celsius"    : round(tm, 1),
        "hairpin_dg"    : hp,
        "homodimer_dg"  : hd,
        "gc_ok"         : gc_ok,
        "tm_ok"         : tm_ok,
        "no_fold"       : fold_ok,
        "no_dimer"      : dimer_ok,
        "PASS"          : pass_all,
    }

# ══════════════════════════════════════════════════════════════════════════════
# FETCH NCBI
# ══════════════════════════════════════════════════════════════════════════════
def safe_fetch(query, retmax=MAX_SEQS):
    for attempt in range(3):
        try:
            time.sleep(0.5)
            h = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
            r = Entrez.read(h); h.close()
            ids = r.get("IdList", [])
            if not ids:
                return None
            time.sleep(0.5)
            h    = Entrez.efetch(db="nucleotide", id=",".join(ids),
                                 rettype="fasta", retmode="text")
            data = h.read(); h.close()
            return data
        except Exception as e:
            print(f"    ⚠ Tentativa {attempt+1}/3: {e}")
            time.sleep(2 * (attempt + 1))
    return None

def fetch_reference(accession: str) -> list:
    """Descarrega a sequência de referência RefSeq pelo accession."""
    try:
        time.sleep(0.5)
        h    = Entrez.efetch(db="nucleotide", id=accession,
                             rettype="fasta", retmode="text")
        data = h.read(); h.close()
        recs = list(SeqIO.parse(StringIO(data), "fasta"))
        return recs
    except Exception as e:
        print(f"    ⚠ Não foi possível descarregar referência {accession}: {e}")
        return []

def filter_sequences(fasta_text, gene_name, min_len, max_len):
    lines = fasta_text.splitlines()
    start = next((i for i, l in enumerate(lines)
                  if l.startswith((">"," #","!",";"))), 0)
    cleaned = "\n".join(lines[start:]).strip()
    if not cleaned:
        return []
    records, seen = [], set()
    for rec in SeqIO.parse(StringIO(cleaned), "fasta"):
        if rec.id in seen:
            continue
        if not (min_len <= len(rec.seq) <= max_len):
            continue
        seen.add(rec.id)
        rec.id          = f"{rec.id.split('.')[0]}_{gene_name}"
        rec.description = ""
        records.append(rec)
    return records

# ══════════════════════════════════════════════════════════════════════════════
# MAFFT
# ══════════════════════════════════════════════════════════════════════════════
def run_mafft(input_fasta: Path, output_fasta: Path, n_seqs: int):
    algo   = "--auto" if n_seqs <= 200 else "--retree 1"
    cmd    = [MAFFT_CMD, algo, "--thread", "4", "--preservecase", str(input_fasta)]
    print(f"    ▶ mafft {algo}  ({n_seqs} seqs)")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            print(f"    ✗ MAFFT erro:\n{r.stderr[:300]}")
            return False
        output_fasta.write_text(r.stdout, encoding="utf-8")
        return True
    except FileNotFoundError:
        print("✗ MAFFT não encontrado no path definido")
        return False
    except subprocess.TimeoutExpired:
        print("    ✗ MAFFT timeout (300s)."); return False

# ══════════════════════════════════════════════════════════════════════════════
# ANÁLISE DE CONSERVAÇÃO
# ══════════════════════════════════════════════════════════════════════════════
def analyse_conservation(aligned_fasta: Path, gene_name: str, target_info: dict):
    records = list(SeqIO.parse(aligned_fasta, "fasta"))
    if len(records) < 1:
        return [], []

    n_seqs  = len(records)
    aln_len = len(records[0].seq)
    cons_thr = target_info.get("cons_threshold", 0.85)
    gap_thr  = target_info.get("gap_threshold",  0.20)

    # Conservação por posição
    conservation = []
    for i in range(aln_len):
        col       = [str(r.seq[i]).upper() for r in records]
        top_count = Counter(col).most_common(1)[0][1]
        gap_freq  = col.count("-") / n_seqs
        conservation.append((top_count / n_seqs, gap_freq))

    # Sliding window para janelas candidatas
    all_windows = []
    for start in range(aln_len):
        for length in range(PROBE_LEN_MIN, PROBE_LEN_MAX + 1):
            end = start + length
            if end > aln_len:
                break
            w        = conservation[start:end]
            avg_cons = sum(c for c, g in w) / length
            avg_gap  = sum(g for c, g in w) / length
            if avg_cons < cons_thr or avg_gap > gap_thr:
                continue
            # Sequência consenso (sem gaps)
            consensus = ""
            for i in range(start, end):
                col = [str(r.seq[i]).upper() for r in records
                       if str(r.seq[i]) != "-"]
                consensus += Counter(col).most_common(1)[0][0] if col else "N"
            if "N" * 3 in consensus:       # rejeitar se ≥3 N consecutivos
                continue
            all_windows.append({
                "start"            : start,
                "end"              : end,
                "length"           : length,
                "avg_conservation" : round(avg_cons, 3),
                "avg_gap"          : round(avg_gap, 3),
                "consensus"        : consensus,
            })

    # Remover sobreposições (manter melhor)
    all_windows.sort(key=lambda x: -x["avg_conservation"])
    best, used = [], set()
    for w in all_windows:
        if len(set(range(w["start"], w["end"])) & used) < w["length"] // 2:
            best.append(w)
            used.update(range(w["start"], w["end"]))
        if len(best) >= 15:
            break
    return conservation, best

# ══════════════════════════════════════════════════════════════════════════════
# PROBES A PARTIR DE SEQUÊNCIA ÚNICA (fallback rmpM)
# ══════════════════════════════════════════════════════════════════════════════
def probes_from_single_seq(seq: str, gene_name: str) -> list:
    """
    Quando não há alinhamento disponível (ex: rmpM com <2 seqs),
    gera probes por sliding window sobre a sequência de referência.
    Conservação = 1.0 (trivialmente, só há 1 sequência).
    """
    seq = seq.upper().replace("-", "")
    windows = []
    for start in range(len(seq)):
        for length in range(PROBE_LEN_MIN, PROBE_LEN_MAX + 1):
            end = start + length
            if end > len(seq):
                break
            candidate = seq[start:end]
            if "N" * 3 in candidate:
                continue
            gc = (candidate.count("G") + candidate.count("C")) / length
            windows.append({
                "start"            : start,
                "end"              : end,
                "length"           : length,
                "avg_conservation" : 1.0,    # só 1 seq — conservação trivial
                "avg_gap"          : 0.0,
                "consensus"        : candidate,
            })
    # Filtrar por GC e remover sobreposições
    windows = [w for w in windows if GC_MIN <= (w["consensus"].count("G")+w["consensus"].count("C"))/w["length"] <= GC_MAX]
    windows.sort(key=lambda x: x["start"])
    best, used = [], set()
    for w in windows:
        if len(set(range(w["start"], w["end"])) & used) < w["length"] // 2:
            best.append(w)
            used.update(range(w["start"], w["end"]))
        if len(best) >= 15:
            break
    return best

# ══════════════════════════════════════════════════════════════════════════════
# SCORING DAS JANELAS (adicionar métricas ViruScope)
# ══════════════════════════════════════════════════════════════════════════════
def score_windows(windows: list) -> list:
    """Adiciona scores ViruScope a cada janela candidata."""
    scored = []
    for w in windows:
        metrics = score_probe(w["consensus"])
        w.update(metrics)
        scored.append(w)
    # Ordenar: PASS primeiro, depois por conservação
    scored.sort(key=lambda x: (-int(x.get("PASS", False)),
                                -x["avg_conservation"],
                                -x.get("tm_celsius", 0)))
    return scored

# ══════════════════════════════════════════════════════════════════════════════
# OUTPUTS
# ══════════════════════════════════════════════════════════════════════════════
def write_outputs(gene_dir: Path, gene_name: str, target_info: dict,
                  n_seqs: int, conservation: list, windows: list,
                  single_seq_mode: bool = False):

    cons_thr = target_info.get("cons_threshold", 0.85)
    gap_thr  = target_info.get("gap_threshold",  0.20)

    # ── 1. conserved_report.txt ───────────────────────────────────────────
    rpt = gene_dir / "conserved_report.txt"
    with open(rpt, "w", encoding="utf-8") as f:
        f.write(f"{'═'*70}\n")
        f.write(f"  RELATÓRIO: {gene_name.upper()}  |  Grupo {target_info['group']}\n")
        f.write(f"  {target_info['organism']}\n")
        f.write(f"  {target_info['notes']}\n")
        f.write(f"{'═'*70}\n\n")
        if single_seq_mode:
            f.write(f"  ⚠ MODO SEQUÊNCIA ÚNICA — poucos homólogos disponíveis no NCBI.\n")
            f.write(f"  Probes geradas por sliding window sobre a sequência de referência.\n")
            f.write(f"  Conservação = 1.0 (trivial). BLAST de especificidade é OBRIGATÓRIO.\n\n")
        f.write(f"Sequências usadas      : {n_seqs}\n")
        if conservation:
            avg_g = sum(c for c,g in conservation) / len(conservation)
            h90   = sum(1 for c,g in conservation if c >= 0.90)
            f.write(f"Comprimento alinhamento: {len(conservation)} posições\n")
            f.write(f"Conservação média      : {avg_g:.3f}\n")
            f.write(f"Posições ≥ 90% cons    : {h90}/{len(conservation)} ({100*h90/len(conservation):.1f}%)\n")
        f.write(f"\n{'─'*70}\n")
        f.write(f"  PROBES CANDIDATAS — ViruScope Scoring (primer3-py)\n")
        f.write(f"  Filtro alinhamento: cons ≥ {cons_thr:.2f}, gap ≤ {gap_thr:.0%}\n")
        f.write(f"  Filtro ViruScope  : Tm {TM_MIN}–{TM_MAX}°C | GC {GC_MIN:.0%}–{GC_MAX:.0%} | hairpin ΔG > {HAIRPIN_DG_MAX} | dimer ΔG > {HOMODIMER_DG_MAX}\n")
        f.write(f"  Critério no-fold  : hairpin ΔG > {HAIRPIN_DG_MAX} kcal/mol (ssDNA para hibridação directa)\n")
        f.write(f"{'─'*70}\n\n")

        pass_count = sum(1 for w in windows if w.get("PASS", False))
        f.write(f"Total de candidatas: {len(windows)}   PASS todas as métricas: {pass_count}\n\n")

        for i, w in enumerate(windows, 1):
            status = "✔ PASS" if w.get("PASS") else "✘ FAIL"
            flags  = []
            if not w.get("gc_ok"):    flags.append(f"GC={w['gc_content']:.0%} fora 40-60%")
            if not w.get("tm_ok"):    flags.append(f"Tm={w['tm_celsius']:.1f}°C fora {TM_MIN}-{TM_MAX}")
            if not w.get("no_fold"):  flags.append(f"hairpin ΔG={w['hairpin_dg']} ≤ {HAIRPIN_DG_MAX} (faz fold!)")
            if not w.get("no_dimer"): flags.append(f"dimer ΔG={w['homodimer_dg']} ≤ {HOMODIMER_DG_MAX}")
            flag_str = "  ← " + "; ".join(flags) if flags else ""

            f.write(f"  #{i:02d} [{status}]  5'-{w['consensus']}-3'\n")
            f.write(f"       pos {w['start']:>5}–{w['end']:>5}  len={w['length']} nt"
                    f"  cons={w['avg_conservation']:.3f}  gap={w['avg_gap']:.3f}\n")
            f.write(f"       Tm={w.get('tm_celsius','N/A'):.1f}°C  "
                    f"GC={w['gc_content']:.0%}  "
                    f"hairpin ΔG={w.get('hairpin_dg','N/A')} kcal/mol  "
                    f"dimer ΔG={w.get('homodimer_dg','N/A')} kcal/mol"
                    f"{flag_str}\n\n")

        f.write(f"{'─'*70}\n  BLAST (especificidade obrigatória)\n{'─'*70}\n")
        f.write(f"  URL: https://blast.ncbi.nlm.nih.gov/Blast.cgi\n")
        f.write(f"  Excluir: {', '.join(target_info.get('cross_react', []))}\n\n")
        f.write(f"{'─'*70}\n  PRÓXIMOS PASSOS\n{'─'*70}\n")
        f.write(f"  1. Abrir aligned.fasta em Geneious/AliView\n")
        f.write(f"  2. Importar {gene_name}_probes_scored.tsv no ViruScope (parse in-silico primers)\n")
        f.write(f"  3. NUPACK: ΔG < -12 kcal/mol, ensemble defect < 0.10\n")
        f.write(f"  4. BLAST especificidade\n")
        f.write(f"  5. Probe PASS final → Fase 1 workflow (3dDNA → docking)\n")

    # ── 2. probe_windows_scored.tsv ───────────────────────────────────────
    tsv = gene_dir / f"{gene_name}_probes_scored.tsv"
    cols = ["rank", "gene", "start", "end", "length",
            "avg_conservation", "avg_gap", "gc_content",
            "tm_celsius", "hairpin_dg", "homodimer_dg",
            "gc_ok", "tm_ok", "no_fold", "no_dimer", "PASS",
            "consensus_5to3"]
    with open(tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for i, w in enumerate(windows, 1):
            row = [
                i, gene_name, w["start"], w["end"], w["length"],
                w["avg_conservation"], w["avg_gap"], w["gc_content"],
                w.get("tm_celsius", "N/A"), w.get("hairpin_dg", "N/A"),
                w.get("homodimer_dg", "N/A"),
                w.get("gc_ok", ""), w.get("tm_ok", ""),
                w.get("no_fold", ""), w.get("no_dimer", ""),
                w.get("PASS", ""),
                w["consensus"],
            ]
            f.write("\t".join(str(x) for x in row) + "\n")

    # ── 3. ViruScope-ready FASTA (só PASS, ou top-5 se nenhum passar) ────
    pass_windows = [w for w in windows if w.get("PASS")]
    to_export    = pass_windows if pass_windows else windows[:5]
    vs_fasta = gene_dir / f"{gene_name}_viroscope_probes.fasta"
    with open(vs_fasta, "w", encoding="utf-8") as f:
        for i, w in enumerate(to_export, 1):
            status = "PASS" if w.get("PASS") else "FAIL"
            f.write(f">{gene_name}_probe_{i:02d}_{status}_"
                    f"pos{w['start']}-{w['end']}_"
                    f"Tm{w.get('tm_celsius','NA')}_"
                    f"GC{w['gc_content']:.2f}_"
                    f"hp{w.get('hairpin_dg','NA')}\n")
            f.write(f"{w['consensus']}\n")

    return rpt, tsv, vs_fasta

# ══════════════════════════════════════════════════════════════════════════════
# RESUMO FINAL CONSOLIDADO (todos os genes)
# ══════════════════════════════════════════════════════════════════════════════
def write_final_summary(all_results: dict):
    """
    Escreve um TSV único com TOP 5 probes por gene.
    Se não houver PASS, usa as melhores disponíveis com nota.
    """
    summary_path = OUTPUT_DIR / "FINAL_PROBES_SUMMARY.tsv"
    rows = []

    for gene, data in all_results.items():

        if not data or not data.get("windows"):
            rows.append({
                "gene": gene,
                "organism": TARGETS[gene]["organism"],
                "group": TARGETS[gene]["group"],
                "status": "FALHOU",
                "probe_5to3": "",
                "length": "",
                "tm_celsius": "",
                "gc_content": "",
                "hairpin_dg": "",
                "avg_conservation": "",
                "PASS": "",
                "notes": data.get("error", ""),
            })
            continue

        windows = data["windows"]

        # 🔥 TOP 5 (já ordenado no pipeline)
        best_list = windows[:5]

        # verifica se alguma PASS existe no top 5
        any_pass = any(w.get("PASS") for w in best_list)

        note = "" if any_pass else "⚠ nenhuma PASS — top 5 usado"

        for w in best_list:
            rows.append({
                "gene": gene,
                "organism": TARGETS[gene]["organism"],
                "group": TARGETS[gene]["group"],
                "status": "PASS" if w.get("PASS") else "BEST_AVAILABLE",
                "probe_5to3": w["consensus"],
                "length": w["length"],
                "tm_celsius": w.get("tm_celsius", ""),
                "gc_content": w.get("gc_content", ""),
                "hairpin_dg": w.get("hairpin_dg", ""),
                "avg_conservation": w["avg_conservation"],
                "PASS": w.get("PASS", ""),
                "notes": note,
            })

    cols = [
        "gene","organism","group","status","probe_5to3","length",
        "tm_celsius","gc_content","hairpin_dg","avg_conservation","PASS","notes"
    ]

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    print(f"\n  {'─'*60}")
    print(f"  LISTA FINAL DE PROBES (TOP 5 por gene)")
    print(f"  {'─'*60}")

    for r in rows:
        flag = "✔ PASS" if r["status"] == "PASS" else "⚠ BEST"

        tm = r.get("tm_celsius")
        tm_str = f"{tm:.1f}" if isinstance(tm, (int, float)) else "NA"

        gc = r.get("gc_content")
        gc_str = f"{gc:.2f}" if isinstance(gc, (int, float)) else "NA"

        hp = r.get("hairpin_dg")
        hp_str = f"{hp:.2f}" if isinstance(hp, (int, float)) else "NA"

        print(
            f"  {flag}  {r['gene']:6s}  "
            f"Tm={tm_str:>5}°C  "
            f"GC={gc_str:>5}  "
            f"hp={hp_str:>6}  "
            f"5'-{str(r['probe_5to3'])[:30]}...-3'"
        )

    print(f"\n  Resumo guardado em: {summary_path}")
    return summary_path

# ══════════════════════════════════════════════════════════════════════════════
# PIPELINE PRINCIPAL POR TARGET
# ══════════════════════════════════════════════════════════════════════════════
def process_target(gene_name: str, target_info: dict, skip_fetch: bool = False):
    print(f"\n{'━'*60}")
    print(f"  {gene_name.upper()}  |  {target_info['organism']}  |  Grupo {target_info['group']}")
    print(f"{'━'*60}")

    gene_dir      = OUTPUT_DIR / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    raw_fasta     = gene_dir / "raw_sequences.fasta"
    aligned_fasta = gene_dir / "aligned.fasta"
    single_seq_mode = False

    min_len = target_info.get("min_len", 200)
    max_len = target_info.get("max_len", 5000)

    # ── 1. FETCH ──────────────────────────────────────────────────────────
    print("  [1/5] Descarregando sequências do NCBI...")
    if skip_fetch and raw_fasta.exists():
        print(f"    ℹ --skip-fetch: a usar {raw_fasta}")
        records = list(SeqIO.parse(raw_fasta, "fasta"))
    else:
        data = safe_fetch(target_info["query_main"])
        if not data:
            print("    ⚠ Query principal sem resultados — a tentar alternativa...")
            data = safe_fetch(target_info["query_alt"])
        if not data and "query_fallback" in target_info:
            print("    ⚠ Query alternativa sem resultados — a tentar fallback...")
            data = safe_fetch(target_info["query_fallback"])

        if data:
            records = filter_sequences(data, gene_name, min_len, max_len)
        else:
            records = []

        if len(records) < 2:
            if target_info.get("allow_single_seq") and len(records) >= 1:
                print(f"    ⚠ Apenas {len(records)} seq encontrada — modo sequência única (rmpM).")
                single_seq_mode = True
            elif len(records) < 1:
                # Tentar descarregar referência
                print(f"    ⚠ Sem resultados — a tentar descarregar referência {target_info['ref_accession']}...")
                ref_recs = fetch_reference(target_info["ref_accession"])
                if ref_recs:
                    records = ref_recs
                    single_seq_mode = True
                    print(f"    ✔ Referência descarregada ({len(str(ref_recs[0].seq))} bp).")
                else:
                    print(f"    ✗ Sem dados para {gene_name}. Descarregar manualmente:")
                    print(f"    → https://www.ncbi.nlm.nih.gov/nucleotide")
                    print(f"    → Guardar como {raw_fasta} e re-correr com --skip-fetch")
                    return {"ok": False, "error": "sem dados NCBI"}
            else:
                print(f"    ✗ Apenas {len(records)} sequência — insuficiente para alinhamento.")
                print(f"    → Adicionar mais sequências manualmente a {raw_fasta} e usar --skip-fetch")
                return {"ok": False, "error": f"só {len(records)} seq após filtragem {min_len}-{max_len}bp"}

        if records:
            SeqIO.write(records, raw_fasta, "fasta")

    print(f"    ✔ {len(records)} sequências  |  modo: {'sequência única' if single_seq_mode else 'alinhamento múltiplo'}")

    # ── 2. ALINHAR ────────────────────────────────────────────────────────
    conservation = []
    if single_seq_mode:
        print("  [2/5] Alinhamento ignorado (sequência única).")
        seq_str = str(records[0].seq)
    else:
        print(f"  [2/5] Alinhamento MAFFT ({len(records)} sequências)...")
        if not run_mafft(raw_fasta, aligned_fasta, len(records)):
            return {"ok": False, "error": "MAFFT falhou"}
        aln_records = list(SeqIO.parse(aligned_fasta, "fasta"))
        if not aln_records:
            return {"ok": False, "error": "MAFFT output vazio"}
        print(f"    ✔ {len(aln_records)} seqs × {len(aln_records[0].seq)} posições")
        seq_str = None   # não necessário em modo múltiplo

    # ── 3. CONSERVAÇÃO + JANELAS ──────────────────────────────────────────
    print("  [3/5] Análise de conservação e janelas candidatas...")
    if single_seq_mode:
        windows = probes_from_single_seq(seq_str, gene_name)
        print(f"    ✔ {len(windows)} janelas da referência (sliding window, GC {GC_MIN:.0%}–{GC_MAX:.0%})")
    else:
        conservation, windows = analyse_conservation(aligned_fasta, gene_name, target_info)
        cons_thr = target_info.get("cons_threshold", 0.85)
        gap_thr  = target_info.get("gap_threshold",  0.20)
        print(f"    ✔ {len(windows)} janelas candidatas (cons ≥ {cons_thr:.2f}, gap ≤ {gap_thr:.0%})")

    if not windows:
        print("    ✗ Nenhuma janela encontrada com os critérios actuais.")
        print(f"    → Tentar reduzir cons_threshold ({target_info.get('cons_threshold',0.85)}) em TARGETS")
        return {"ok": False, "error": "sem janelas"}

    # ── 4. SCORING VIROSCOPE ──────────────────────────────────────────────
    print("  [4/5] Scoring ViruScope (Tm, GC, hairpin ΔG, dimer ΔG)...")
    windows = score_windows(windows)
    pass_n  = sum(1 for w in windows if w.get("PASS"))
    print(f"    ✔ {pass_n}/{len(windows)} probes passam todos os critérios")
    if pass_n == 0:
        print(f"    ⚠ Nenhuma probe PASS — verificar thresholds ou ajustar comprimento")

    # ── 5. OUTPUTS ────────────────────────────────────────────────────────
    print("  [5/5] A escrever outputs...")
    rpt, tsv, vs = write_outputs(
        gene_dir, gene_name, target_info,
        len(records), conservation, windows, single_seq_mode)
    print(f"    ✔ {rpt.name}")
    print(f"    ✔ {tsv.name}  ← abrir no Excel / ViruScope")
    print(f"    ✔ {vs.name}   ← importar no ViruScope (parse in-silico primers)")
    if not single_seq_mode:
        print(f"    ✔ aligned.fasta  ← Geneious / AliView")

    return {"ok": True, "windows": windows}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    skip_fetch = "--skip-fetch" in sys.argv
    args       = [a for a in sys.argv[1:] if not a.startswith("--")]

    print(f"  ✔ MAFFT: {MAFFT_CMD}")
    print(f"  ✔ primer3-py: {'sim' if HAS_PRIMER3 else 'NÃO (scoring básico)'}")
    if skip_fetch:
        print("  ℹ --skip-fetch activo")

    targets_to_run = args if args else list(TARGETS.keys())
    invalid = [t for t in targets_to_run if t not in TARGETS]
    if invalid:
        print(f"✗ Targets desconhecidos: {invalid}")
        print(f"  Disponíveis: {list(TARGETS.keys())}")
        sys.exit(1)

    OUTPUT_DIR.mkdir(exist_ok=True)
    print(f"\n{'═'*60}")
    print(f"  GFET Probe Pipeline v3 — Fase 0 + ViruScope Scoring")
    print(f"  Targets: {targets_to_run}")
    print(f"{'═'*60}")

    all_results = {}
    for gene in targets_to_run:
        all_results[gene] = process_target(gene, TARGETS[gene], skip_fetch)

    # Resumo por target
    print(f"\n{'═'*60}")
    print("  RESUMO")
    print(f"{'═'*60}")
    for gene, res in all_results.items():
        ok = res.get("ok") if res else False
        print(f"  {'✔' if ok else '✗'}  {gene:10s}  {TARGETS[gene]['organism']}")

    # Lista final consolidada
    write_final_summary(all_results)

if __name__ == "__main__":
    main()
