"""
GFET Probe Design Pipeline
Etapas: NCBI → MAFFT → conservação → primer3 → seqfold → export Colab (Boltz-2)

Referências:
  SantaLucia & Hicks 2004  — Tm nearest-neighbour, parâmetros primer3
  Wetmur 1991              — comprimento da probe (18–28 nt)
  Zadeh et al. 2011        — limiar seqfold (estrutura secundária)
  IDT OligoAnalyzer        — GC 40–65%, homodimer ΔG > -5 kcal/mol
  Stover et al. 2000       — P. aeruginosa GC genómico ~67% (PAO1)
─────────────────────────────────────────────────────────────────────────────
"""

import csv, itertools, re, shutil, subprocess, time, yaml, zipfile
from dataclasses import dataclass, asdict
from io import StringIO
from pathlib import Path
from typing import Optional

import pandas as pd
import primer3
from Bio import AlignIO, Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "pg45861@uminho.pt"

# ── Parâmetros primer3 (SantaLucia & Hicks 2004) ────────────────────────────
P3_MV_CONC  = 50.0    # [Na+] mM
P3_DV_CONC  = 0.0     # [Mg2+] mM
P3_DNTP     = 0.0     # [dNTP] mM
P3_DNA_CONC = 250.0   # [oligo] nM
P3_TEMP_C   = 37.0    # temperatura de referência

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR     = Path(__file__).parent
ALIGN_DIR    = BASE_DIR / "output" / "alignments"
COLAB_DIR    = BASE_DIR / "output" / "colab_boltz2"
BEATRIZ_XLSX = BASE_DIR / "data" / "BeatrizMasterThesis_ProbesData_IPLEXMED_ENERGIAS_CAPS.xlsx"

_MAFFT_CANDIDATES = [
    str(BASE_DIR / "MAFFT" / "mafft-win" / "mafft.bat"),
    str(BASE_DIR / "MAFFT" / "mafft-7.526-win64-signed" / "mafft-win" / "mafft.bat"),
    shutil.which("mafft"),
    shutil.which("mafft.bat"),
]
MAFFT_BIN = Path(next(p for p in _MAFFT_CANDIDATES if p and Path(p).exists()))

for _d in [ALIGN_DIR, COLAB_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# ── Targets ──────────────────────────────────────────────────────────────────
# Chaves opcionais por gene sobrepõem os DEFAULTS:
#   tm_min, tm_max, gc_min, gc_max, hp_min, dimer_min,
#   cons_min, gap_max, probe_len_min, probe_len_max, seqfold_dg_max
TARGETS = {
    "nuc": {
        "organism": "Staphylococcus aureus", "group": "A",
        "ncbi_query":     'nuc[Gene Name] AND "Staphylococcus aureus"[Organism] AND 200:3000[Sequence Length]',
        "ncbi_query_alt": '"Staphylococcus aureus"[Organism] AND nuc[Title] AND 200:2000[Sequence Length]',
        "refseq_fallback": "NC_002951",
        "min_len": 200, "max_len": 3000,
        "tm_min": 52.0,          # gene AT-rico — SantaLucia & Hicks 2004
        "gc_min": 0.38,          # AT-rico — IDT OligoAnalyzer
        "gc_max": 0.60,
    },
    "rmpM": {
        "organism": "Neisseria meningitidis", "group": "A",
        "ncbi_query":        'rmpM[Gene Name] AND "Neisseria meningitidis"[Organism] AND 100:1200[Sequence Length]',
        "ncbi_query_alt":    '("Neisseria meningitidis"[Organism]) AND (rmpM[Title] OR "class 4 outer membrane protein"[Title]) AND 100:1200[Sequence Length]',
        "ncbi_query_fallback": '"Neisseria meningitidis"[Organism] AND "outer membrane protein"[Title] AND 100:1200[Sequence Length]',
        "refseq_fallback": "NC_003112",
        "min_len": 100, "max_len": 1200,
        "gc_max": 0.65,           # IDT OligoAnalyzer guidelines
        "cons_min": 0.80,         # poucos homólogos NCBI
        "allow_single_seq": True,
    },
    "lytA": {
        "organism": "Streptococcus pneumoniae", "group": "B",
        "ncbi_query":     'lytA[Gene Name] AND "Streptococcus pneumoniae"[Organism] AND 700:1300[Sequence Length]',
        "ncbi_query_alt": '"Streptococcus pneumoniae"[Organism] AND lytA[Title] AND 700:1300[Sequence Length]',
        "refseq_fallback": "NC_003028",
        "min_len": 700, "max_len": 1300,
        "cons_min": 0.70,         # diversidade alélica — Whatmore et al. 2000
        "gap_max":  0.40,
    },
    "oprL": {
        "organism": "Pseudomonas aeruginosa", "group": "B",
        "ncbi_query":     'oprL[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 300:2000[Sequence Length]',
        "ncbi_query_alt": '"Pseudomonas aeruginosa"[Organism] AND oprL[Title] AND 300:2000[Sequence Length]',
        "refseq_fallback": "NC_002516",
        "min_len": 300, "max_len": 2000,
        "gc_max": 0.70,           # P. aeruginosa ~67% GC — Stover et al. 2000
    },
    "algD": {
        "organism": "Pseudomonas aeruginosa", "group": "B",
        "ncbi_query":     'algD[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 500:2500[Sequence Length]',
        "ncbi_query_alt": '"Pseudomonas aeruginosa"[Organism] AND algD[Title] AND 500:2500[Sequence Length]',
        "refseq_fallback": "NC_002516",
        "min_len": 500, "max_len": 2500,
        "gc_max": 0.70,           # P. aeruginosa ~67% GC — Stover et al. 2000
    },
    "frdB": {
        "organism": "Haemophilus influenzae", "group": "B",
        "ncbi_query":     'frdB[Gene Name] AND "Haemophilus influenzae"[Organism] AND 200:2000[Sequence Length]',
        "ncbi_query_alt": '"Haemophilus influenzae"[Organism] AND frdB[Title] AND 200:2000[Sequence Length]',
        "refseq_fallback": "NC_000907",
        "min_len": 200, "max_len": 2000,
        "gc_min": 0.38,           # H. influenzae AT-rico — IDT OligoAnalyzer
        "gc_max": 0.60,
    },
}

# ── Limiares globais ─────────────────────────────────────────────────────────
DEFAULTS = {
    "tm_min":          53.0,   # °C — SantaLucia & Hicks 2004
    "tm_max":          72.0,
    "gc_min":          0.40,   # IDT OligoAnalyzer
    "gc_max":          0.65,
    "hp_min":         -2.0,    # kcal/mol — hairpin ΔG mínimo
    "dimer_min":      -5.0,    # kcal/mol — homodimer ΔG mínimo
    "cons_min":        0.85,   # conservação mínima no alinhamento
    "gap_max":         0.20,   # gap máximo no alinhamento
    "probe_len_min":   18,     # nt — Wetmur 1991
    "probe_len_max":   28,
    "window_step":      3,
    "max_seqs":        20,     # sequências a descarregar do NCBI
    "min_len":        100,     # bp — comprimento mínimo da sequência alvo
    "max_len":       5000,
    "seqfold_dg_max": -6.0,   # kcal/mol — ΔG MFE máximo aceite (Zadeh et al. 2011)
}

def cfg(gene_key: str, param: str):
    return TARGETS[gene_key].get(param, DEFAULTS[param])

# ── Dataclass ────────────────────────────────────────────────────────────────
@dataclass
class Probe:
    gene:               str
    organism:           str
    group:              str
    probe_id:           str   = ""
    sequence:           str   = ""
    pos_start:          int   = 0
    pos_end:            int   = 0
    tm:                 float = 0.0
    gc:                 float = 0.0
    hairpin_dg:         Optional[float] = None
    homodimer_dg:       Optional[float] = None
    seqfold_dg:         Optional[float] = None
    seqfold_paired_frac: Optional[float] = None
    pass_basic:         bool  = False
    pass_seqfold:       bool  = False
    fonte:              str   = "romeu"
    notes:              str   = ""

_ORG_ABBR = {
    "Staphylococcus aureus":    "Saur",
    "Neisseria meningitidis":   "Nmen",
    "Streptococcus pneumoniae": "Spne",
    "Pseudomonas aeruginosa":   "Paer",
    "Haemophilus influenzae":   "Hinf",
}

def build_probe_id(probe: Probe, index: int) -> str:
    org  = _ORG_ABBR.get(probe.organism, probe.organism[:4].capitalize())
    hp   = f"{probe.hairpin_dg:.2f}" if probe.hairpin_dg is not None else "NA"
    stat = "PASS" if probe.pass_basic else "FAIL"
    return (f"p{index:03d}_{org}_{probe.gene}"
            f"_pos{probe.pos_start}-{probe.pos_end}"
            f"_Tm{probe.tm:.1f}_GC{int(round(probe.gc*100))}_hp{hp}_{stat}")

# ── 1. NCBI ──────────────────────────────────────────────────────────────────
def fetch_sequences(gene_key: str) -> list[SeqRecord]:
    t     = TARGETS[gene_key]
    max_s = DEFAULTS["max_seqs"]
    min_l = t.get("min_len", DEFAULTS["min_len"])
    max_l = t.get("max_len", DEFAULTS["max_len"])
    print(f"  [1/5] Descarregando sequências do NCBI...")

    def _run_query(query: str) -> list[SeqRecord]:
        for attempt in range(3):
            try:
                time.sleep(0.5)
                h   = Entrez.esearch(db="nucleotide", term=query, retmax=max_s)
                ids = Entrez.read(h)["IdList"]; h.close()
                if not ids:
                    return []
                time.sleep(0.5)
                h   = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
                raw = "\n".join(l for l in h.read().splitlines()
                                if not l.startswith(("!", "#", ";")))
                h.close()
                recs = list(SeqIO.parse(StringIO(raw), "fasta"))
                return [r for r in recs if min_l <= len(r.seq) <= max_l]
            except Exception as e:
                print(f"    ⚠ Tentativa {attempt+1}/3: {e}")
                time.sleep(2 * (attempt + 1))
        return []

    records = _run_query(t["ncbi_query"])

    if len(records) < 2 and "ncbi_query_alt" in t:
        print("    ⚠ Query principal insuficiente — a tentar alternativa...")
        records = _run_query(t["ncbi_query_alt"])

    if len(records) < 2 and "ncbi_query_fallback" in t:
        print("    ⚠ Query alternativa insuficiente — a tentar fallback...")
        records = _run_query(t["ncbi_query_fallback"])

    if len(records) < 2 and t.get("allow_single_seq") and len(records) == 1:
        print(f"    ✔ 1 sequência  |  modo: sequência única")
        return records

    if not records:
        fb = t.get("refseq_fallback")
        if fb:
            print(f"    ⚠ Sem resultados — a usar referência {fb}")
            try:
                time.sleep(0.5)
                h   = Entrez.efetch(db="nucleotide", id=fb, rettype="fasta", retmode="text")
                raw = "\n".join(l for l in h.read().splitlines()
                                if not l.startswith(("!", "#", ";")))
                h.close()
                records = list(SeqIO.parse(StringIO(raw), "fasta"))
            except Exception as e:
                raise RuntimeError(f"Sem dados para {gene_key}: {e}")
        else:
            raise RuntimeError(f"Sem dados para {gene_key}")

    print(f"    ✔ {len(records)} sequências  |  modo: alinhamento múltiplo")
    return records

# ── 2. MAFFT ─────────────────────────────────────────────────────────────────
def align_mafft(records: list[SeqRecord], gene_key: str) -> AlignIO.MultipleSeqAlignment:
    out_dir = ALIGN_DIR / gene_key
    out_dir.mkdir(exist_ok=True)
    fa_in  = out_dir / "input.fasta"
    fa_out = out_dir / "aligned.fasta"
    if len(records) == 1:
        print(f"  [2/5] Alinhamento ignorado (sequência única).")
        SeqIO.write(records, fa_out, "fasta")
        return AlignIO.read(fa_out, "fasta")
    SeqIO.write(records, fa_in, "fasta")
    print(f"  [2/5] Alinhamento MAFFT ({len(records)} sequências)...")
    result = subprocess.run([str(MAFFT_BIN), "--auto", str(fa_in)],
                            capture_output=True, text=True)
    fa_out.write_text(result.stdout)
    aln = AlignIO.read(fa_out, "fasta")
    print(f"    ✔ {len(aln)} seqs × {aln.get_alignment_length()} posições")
    return aln

# ── 3. Janelas candidatas ────────────────────────────────────────────────────
def _conservation(col: list[str]) -> float:
    bases = [b for b in col if b not in ("-", "N", "n")]
    if not bases:
        return 0.0
    return max(bases.count(b) for b in set(bases)) / len(bases)

def _gap_freq(col: list[str]) -> float:
    return col.count("-") / len(col)

def candidate_windows(aln: AlignIO.MultipleSeqAlignment, gene_key: str) -> list[dict]:
    cons_min = cfg(gene_key, "cons_min")
    gap_max  = cfg(gene_key, "gap_max")
    plen_min = cfg(gene_key, "probe_len_min")
    plen_max = cfg(gene_key, "probe_len_max")
    step     = cfg(gene_key, "window_step")
    aln_len  = aln.get_alignment_length()
    cols     = [[str(rec.seq[i]) for rec in aln] for i in range(aln_len)]
    windows  = []
    for start in range(0, aln_len - plen_min + 1, step):
        for plen in range(plen_min, plen_max + 1):
            end = start + plen
            if end > aln_len:
                break
            w_cols = cols[start:end]
            if any(_gap_freq(c) > gap_max for c in w_cols):
                continue
            cons = sum(_conservation(c) for c in w_cols) / plen
            if cons < cons_min:
                continue
            consensus = "".join(
                max(set(c) - {"-"}, key=lambda b: c.count(b), default="N")
                for c in w_cols
            ).upper()
            windows.append({"start": start, "end": end, "cons": cons, "seq": consensus})
    print(f"  [3/5] ✔ {len(windows)} janelas candidatas (cons ≥ {cons_min}, gap ≤ {gap_max*100:.0f}%)")
    return windows

# ── 4. Scoring básico ────────────────────────────────────────────────────────
def score_probe(seq: str) -> dict:
    kw = dict(mv_conc=P3_MV_CONC, dv_conc=P3_DV_CONC,
              dntp_conc=P3_DNTP, dna_conc=P3_DNA_CONC)
    tm = primer3.calc_tm(seq, **kw)
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    try:
        hp_res = primer3.calc_hairpin(seq, **kw, temp_c=P3_TEMP_C)
        hp = round(hp_res.dg / 1000, 2) if hp_res.structure_found else None
        if hp is not None and not (-50.0 < hp < 50.0):
            hp = None
    except Exception:
        hp = None
    try:
        dm_res = primer3.calc_homodimer(seq, **kw, temp_c=P3_TEMP_C)
        dm = round(dm_res.dg / 1000, 2) if dm_res.structure_found else None
        if dm is not None and not (-50.0 < dm < 50.0):
            dm = None
    except Exception:
        dm = None
    return {"tm": round(tm, 1), "gc": round(gc, 3), "hairpin_dg": hp, "homodimer_dg": dm}

def passes_basic(probe: Probe, gene_key: str) -> bool:
    return (
        cfg(gene_key, "tm_min")    <= probe.tm <= cfg(gene_key, "tm_max") and
        cfg(gene_key, "gc_min")    <= probe.gc <= cfg(gene_key, "gc_max") and
        (probe.hairpin_dg   is None or probe.hairpin_dg   >= cfg(gene_key, "hp_min")) and
        (probe.homodimer_dg is None or probe.homodimer_dg >= cfg(gene_key, "dimer_min"))
    )

def print_screening_funnel(probes: list, gene_key: str) -> None:
    total     = len(probes)
    remaining = probes[:]
    tm_min = cfg(gene_key, "tm_min");  tm_max = cfg(gene_key, "tm_max")
    gc_min = cfg(gene_key, "gc_min");  gc_max = cfg(gene_key, "gc_max")
    hp_min = cfg(gene_key, "hp_min");  dm_min = cfg(gene_key, "dimer_min")
    steps = [
        (f"Tm ({tm_min:.0f}–{tm_max:.0f}°C)",
         lambda p: tm_min <= p.tm <= tm_max),
        (f"GC ({gc_min*100:.0f}–{gc_max*100:.0f}%)",
         lambda p: gc_min <= p.gc <= gc_max),
        (f"Hairpin ΔG > {hp_min} kcal/mol",
         lambda p: p.hairpin_dg is None or p.hairpin_dg >= hp_min),
        (f"Homodimer ΔG > {dm_min} kcal/mol",
         lambda p: p.homodimer_dg is None or p.homodimer_dg >= dm_min),
    ]
    print(f"    {'Critério':<30} {'Passam':>6}  {'Rejeitadas':>10}")
    print(f"    {'─'*30} {'─'*6}  {'─'*10}")
    print(f"    {'Total de janelas candidatas':<30} {total:>6}")
    for label, test in steps:
        ok   = [p for p in remaining if test(p)]
        fail = len(remaining) - len(ok)
        print(f"    ↳ {label:<28} {len(ok):>6}  {fail:>10}")
        remaining = ok
    print(f"    {'─'*30} {'─'*6}  {'─'*10}")
    print(f"    {'✔ Passam todos os critérios':<30} {len(remaining):>6}  {total - len(remaining):>10}")

# ── 5. seqfold ───────────────────────────────────────────────────────────────
def run_seqfold_probe(probe: Probe, gene_key: str) -> Probe:
    try:
        import seqfold as sf
    except ImportError:
        probe.notes += "[seqfold não instalado: pip install seqfold] "
        return probe
    seq = probe.sequence
    try:
        dg      = sf.dg(seq, temp=37.0)
        structs = sf.fold(seq, temp=37.0)
        paired  = {idx for s in structs for pair in s.ij for idx in pair if idx >= 0}
        probe.seqfold_dg          = round(dg, 2)
        probe.seqfold_paired_frac = round(len(paired) / len(seq) if seq else 0.0, 3)
    except Exception as e:
        probe.notes += f"[seqfold erro: {e}] "
        return probe
    probe.pass_seqfold = probe.seqfold_dg >= cfg(gene_key, "seqfold_dg_max")
    return probe

# ── 6. Export para Colab (Boltz-2) ───────────────────────────────────────────
def export_colab_inputs(probes: list, top_n: int) -> Path:
    """
    Gera os ficheiros de input para o Boltz-2 Colab:
      - YAML individual por probe (formato nativo Boltz-2)
      - FASTA consolidado
      - boltz2_inputs.zip → upload directo ao Colab batch
    """
    yaml_dir = COLAB_DIR / "yamls"
    yaml_dir.mkdir(exist_ok=True)

    selected: list[Probe] = []
    for gk in TARGETS:
        cands = sorted(
            [p for p in probes if p.gene == gk and p.pass_basic and p.pass_seqfold],
            key=lambda p: p.tm, reverse=True
        )[:top_n]
        selected.extend(cands)

    if not selected:
        print("  ⚠ Nenhuma probe qualificada para export Colab.")
        return COLAB_DIR

    for p in selected:
        (yaml_dir / f"{p.probe_id}.yaml").write_text(
            yaml.dump({"version": 1,
                       "sequences": [{"dna": {"id": ["A"], "sequence": p.sequence}}]}),
            encoding="utf-8"
        )

    with open(COLAB_DIR / "probes_for_boltz2.fasta", "w", encoding="utf-8") as f:
        for p in selected:
            f.write(f">{p.probe_id}\n{p.sequence}\n")

    meta_fields = ["probe_id", "gene", "organism", "sequence",
                   "tm", "gc", "hairpin_dg", "homodimer_dg", "seqfold_dg"]
    with open(COLAB_DIR / "probes_metadata.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=meta_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows([asdict(p) for p in selected])

    zip_path = COLAB_DIR / "boltz2_inputs.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        for yf in sorted(yaml_dir.glob("*.yaml")):
            zf.write(yf, yf.name)

    print(f"\n{'═'*60}")
    print(f"  EXPORT COLAB — Boltz-2")
    print(f"  {len(selected)} probes seleccionadas (top {top_n} por gene, por Tm)")
    print(f"{'═'*60}")
    print(f"  Ficheiros em: {COLAB_DIR}")
    print(f"    boltz2_inputs.zip         ← upload para o Colab batch")
    print(f"    probes_for_boltz2.fasta   ← FASTA alternativo")
    print(f"    probes_metadata.csv       ← metadados (Tm, GC, etc.)")
    print(f"\n  Colab batch (recomendado):")
    print(f"    https://colab.research.google.com/github/AtharvaTilewale/boltz2-notebook/blob/main/batch/Boltz2_Batch_v1_beta.ipynb")
    print(f"\n  Passos:")
    print(f"    1. Abrir o Colab acima")
    print(f"    2. Upload de boltz2_inputs.zip")
    print(f"    3. Correr todas as células")
    print(f"    4. Download do ZIP de resultados (CIF files)")
    print(f"    5. Extrair os CIF para output/structures/")
    print(f"{'═'*60}")
    return zip_path

# ── 7. Outputs por gene ──────────────────────────────────────────────────────
def write_gene_outputs(probes: list[Probe], gene_key: str):
    out_dir = ALIGN_DIR / gene_key
    tsv   = out_dir / f"{gene_key}_probes_scored.tsv"
    fasta = out_dir / f"{gene_key}_viroscope_probes.fasta"
    fields = list(asdict(probes[0]).keys())
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for p in probes:
            w.writerow(asdict(p))
    with open(fasta, "w") as f:
        for p in probes:
            if p.pass_basic:
                f.write(f">{p.probe_id}\n{p.sequence}\n")
    print(f"    ✔ {tsv.name}  ✔ {fasta.name}")

# ── 8. CSV consolidado (Romeu + Beatriz) ─────────────────────────────────────
def _parse_beatriz_xlsx(path: Path) -> list[dict]:
    rows: list[dict] = []
    if not path.exists():
        print(f"  ⚠ Ficheiro Beatriz não encontrado: {path}")
        return rows

    _NORM_ORG = {"Strepptocuccus": "Streptococcus", "Strepptococcus": "Streptococcus"}
    _GROUP_A  = ["Staphylococcus aureus", "Neisseria meningitidis"]

    def _clean_seq(s) -> str:
        if pd.isna(s): return ""
        s = re.sub(r"^[A-Z0-9]+-C\d+-", "", str(s).upper())
        return re.sub(r"[^ATGCUN]", "", s)

    def _parse_float(s) -> Optional[float]:
        if pd.isna(s): return None
        m = re.search(r"[-\d.]+", str(s))
        return float(m.group()) if m else None

    def _parse_gc(s) -> Optional[float]:
        v = _parse_float(s)
        return round(v / 100, 3) if v is not None else None

    def _norm_org(s) -> str:
        if pd.isna(s): return ""
        s = str(s).strip()
        for wrong, right in _NORM_ORG.items():
            s = s.replace(wrong, right)
        return s

    def _extract_gene(s) -> str:
        if pd.isna(s): return ""
        s = str(s).strip()
        return s.split(":")[0].strip() if ":" in s else s.split()[0]

    def _group(org: str) -> str:
        return "A" if any(g in org for g in _GROUP_A) else "B"

    def _thresh(gene: str, param: str):
        gk = gene.lower() if gene.lower() in TARGETS else None
        return TARGETS[gk].get(param, DEFAULTS[param]) if gk else DEFAULTS[param]

    def _pass_basic(tm, gc, hp, gene: str) -> Optional[bool]:
        if tm is None or gc is None: return None
        return bool(
            _thresh(gene, "tm_min") <= tm <= _thresh(gene, "tm_max") and
            _thresh(gene, "gc_min") <= gc <= _thresh(gene, "gc_max") and
            (hp is None or hp >= _thresh(gene, "hp_min"))
        )

    xl = pd.ExcelFile(path)

    _raw1 = pd.read_excel(xl, sheet_name="Nossas sequencias", header=None)
    _hrow = next(
        (i for i, r in _raw1.iterrows()
         if any(str(v).strip().lower() == "sequence" for v in r.values if pd.notna(v))), 2)
    df1 = pd.read_excel(xl, sheet_name="Nossas sequencias", header=_hrow, dtype=str)
    df1.columns = [str(c).strip() for c in df1.columns]
    for _, row in df1.iterrows():
        seq = _clean_seq(row.get("Sequence"))
        if len(seq) < 10: continue
        organism = _norm_org(row.get("species"))
        gene     = _extract_gene(row.get("Gene/STRAIN"))
        tm       = _parse_float(row.get("Tm"))
        gc       = _parse_gc(row.get("GC content"))
        hp       = _parse_float(row.get("Lowest folding energy kcal/mole at 25C"))
        pid_raw  = row.get("Column1")
        pid      = str(pid_raw).strip() if pd.notna(pid_raw) else f"bea{len(rows)+1:03d}"
        rows.append({
            "probe_id": pid, "gene": gene, "organism": organism, "group": _group(organism),
            "sequence": seq, "pos_start": "", "pos_end": "",
            "tm": tm if tm is not None else "", "gc": gc if gc is not None else "",
            "hairpin_dg": hp, "homodimer_dg": "",
            "seqfold_dg": "", "seqfold_paired_frac": "",
            "pass_basic": _pass_basic(tm, gc, hp, gene), "pass_seqfold": "",
            "fonte": "beatriz", "notes": str(row.get("notes", "") or "").strip(),
        })

    _raw2 = pd.read_excel(xl, sheet_name="New probes (papers)", header=None)
    _hrow2 = next(
        (i for i, r in _raw2.iterrows()
         if any("probe" in str(v).lower() or "species" in str(v).lower()
                for v in r.values if pd.notna(v))), 0)
    df2 = pd.read_excel(xl, sheet_name="New probes (papers)", header=_hrow2, dtype=str)
    df2.columns = [str(c).strip() for c in df2.columns]
    for _, row in df2.iterrows():
        seq = re.sub(r"[^ATGCUN]", "",
                     str(row.get("Published probes (5'-3')", "") or "").upper())
        if len(seq) < 10: continue
        organism = _norm_org(row.get("species"))
        gene     = str(row.get("gene", "") or "").strip()
        source   = str(row.get("Source", "") or "").strip()
        notes    = str(row.get("notes", "") or "").strip()
        if source:
            notes = f"{notes} [Source: {source}]".strip("[ ]")
        sc  = score_probe(seq)
        pid = (f"lit{len(rows)+1:03d}_"
               f"{_ORG_ABBR.get(organism, organism[:4].capitalize())}_{gene}")
        rows.append({
            "probe_id": pid, "gene": gene, "organism": organism, "group": _group(organism),
            "sequence": seq, "pos_start": "", "pos_end": "",
            "tm": sc["tm"], "gc": sc["gc"],
            "hairpin_dg": sc["hairpin_dg"], "homodimer_dg": sc["homodimer_dg"],
            "seqfold_dg": "", "seqfold_paired_frac": "",
            "pass_basic": _pass_basic(sc["tm"], sc["gc"], sc["hairpin_dg"], gene),
            "pass_seqfold": "",
            "fonte": "literatura", "notes": notes,
        })

    print(f"  ✔ Beatriz: {len(rows)} probes lidas do XLSX")
    return rows

def write_consolidated_csv(all_probes: list[Probe]):
    out_path     = BASE_DIR / "output" / "FINAL_PROBES_ALL.csv"
    romeu_rows   = [asdict(p) for p in all_probes]
    beatriz_rows = _parse_beatriz_xlsx(BEATRIZ_XLSX)
    all_rows     = romeu_rows + beatriz_rows
    all_rows.sort(key=lambda r: (r["gene"], not r["pass_basic"]))
    fields = list(asdict(Probe("", "", "")).keys())
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"\n  ✔ CSV consolidado: {out_path}  ({len(all_rows)} probes total)")

# ── Pipeline principal ────────────────────────────────────────────────────────
def run_pipeline(run_seqfold: bool = True, colab_top: int = 0):
    print("\n" + "═"*60)
    print("  GFET Probe Pipeline")
    print(f"  Targets: {list(TARGETS.keys())}")
    print("═"*60)

    all_probes: list[Probe] = []
    probe_counter = itertools.count(1)

    for gene_key, t in TARGETS.items():
        print(f"\n{'━'*60}")
        print(f"  {gene_key.upper()}  |  {t['organism']}  |  Grupo {t['group']}")
        print(f"{'━'*60}")

        records = fetch_sequences(gene_key)
        aln     = align_mafft(records, gene_key)
        windows = candidate_windows(aln, gene_key)

        probes: list[Probe] = []
        print(f"  [4/5] Scoring básico (Tm, GC, hairpin, homodimer)...")
        for w in windows:
            sc = score_probe(w["seq"])
            p  = Probe(
                gene=gene_key, organism=t["organism"], group=t["group"],
                sequence=w["seq"], pos_start=w["start"], pos_end=w["end"],
                tm=sc["tm"], gc=sc["gc"],
                hairpin_dg=sc["hairpin_dg"], homodimer_dg=sc["homodimer_dg"],
            )
            p.pass_basic = passes_basic(p, gene_key)
            p.probe_id   = build_probe_id(p, next(probe_counter))
            probes.append(p)

        print_screening_funnel(probes, gene_key)
        n_pass = sum(p.pass_basic for p in probes)

        if run_seqfold:
            print(f"  [5/5] seqfold (ΔG MFE + fracção emparelhada)...")
            for p in probes:
                if p.pass_basic:
                    p = run_seqfold_probe(p, gene_key)
            n_sf = sum(p.pass_seqfold for p in probes if p.pass_basic)
            print(f"    ✔ {n_sf}/{n_pass} passam seqfold  [{n_pass - n_sf} rejeitadas]")

        write_gene_outputs(probes, gene_key)
        all_probes.extend(probes)

    write_consolidated_csv(all_probes)

    print("\n" + "═"*60)
    print("  RESUMO")
    print("═"*60)
    for gk in TARGETS:
        ps = [p for p in all_probes if p.gene == gk]
        nb = sum(p.pass_basic   for p in ps)
        ns = sum(p.pass_seqfold for p in ps)
        print(f"  {gk:<8} basic={nb:<5} seqfold={ns}")

    if colab_top > 0:
        export_colab_inputs(all_probes, top_n=colab_top)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="GFET Probe Pipeline")
    ap.add_argument("--no-seqfold", action="store_true",
                    help="Saltar análise de estrutura secundária (seqfold)")
    ap.add_argument("--colab", type=int, default=0, metavar="N",
                    help="Gerar YAMLs + ZIP para Boltz-2 Colab (top N por gene, ex: --colab 5)")
    args = ap.parse_args()
    run_pipeline(run_seqfold=not args.no_seqfold, colab_top=args.colab)
