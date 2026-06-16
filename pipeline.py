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

import csv, itertools, math, re, shutil, subprocess, time, yaml, zipfile
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
    "tm_min":          53.0,   # °C — SantaLucia & Hicks 2004 (default; ver assay_temp_c)
    "tm_max":          72.0,
    "assay_temp_c":    None,   # se definida (°C), deriva a janela Tm automaticamente
    "tm_offset_min":   15.0,   # Tm_min = T_ensaio + 15  (duplex estável acima da T do ensaio)
    "tm_offset_max":   35.0,   # Tm_max = T_ensaio + 35  (a 37°C → 52–72, ~janela atual)
    "gc_min":          0.40,   # IDT OligoAnalyzer
    "gc_max":          0.60,   # janela ótima 40–60% (PremierBiosoft/UNLV); override 0.70 P. aeruginosa
    "hp_min":         -2.0,    # kcal/mol — hairpin ΔG mínimo
    "dimer_min":      -5.0,    # kcal/mol — homodimer ΔG mínimo
    "cons_min":        0.85,   # conservação mínima no alinhamento
    "gap_max":         0.20,   # gap máximo no alinhamento
    "probe_len_min":   18,     # nt — Wetmur 1991
    "probe_len_max":   28,
    "window_step":      3,
    "max_seqs":        20,     # sequências a descarregar do NCBI
    "min_len":        100,     # bp — pré-filtro grosseiro (sanidade) do comprimento
    "max_len":       5000,
    "len_cluster":     True,   # selecionar automaticamente o cluster de comprimento
    "len_cluster_tol": 0.25,   # ±25% em torno do comprimento dominante
    "seqfold_dg_max": -2.6,   # kcal/mol — P5 da distribuição observada (fixo/reprodutível); era -6.0
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
    conservation:       Optional[float] = None   # PPI (ViruScope), fração 0–1
    hairpin_dg:         Optional[float] = None
    homodimer_dg:       Optional[float] = None
    seqfold_dg:         Optional[float] = None
    seqfold_paired_frac: Optional[float] = None
    nofold_score:       Optional[float] = None   # No-fold score (ViruScope), 0–100
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
def _parse_efetch_fasta(raw: str) -> list[SeqRecord]:
    """Parse robusto do FASTA do efetch, independente da versão do Biopython:
    descarta tudo antes do primeiro '>' (o efetch devolve frequentemente uma linha
    em branco inicial que o parser 'fasta' estrito rejeitaria) e usa o formato
    'fasta' universal — evita 'fasta-pearson', que não existe em Biopython < 1.85."""
    i = raw.find(">")
    if i < 0:
        return []
    return list(SeqIO.parse(StringIO(raw[i:]), "fasta"))

def fetch_sequences(gene_key: str) -> list[SeqRecord]:
    t     = TARGETS[gene_key]
    max_s = cfg(gene_key, "max_seqs")          # configurável (CLI/por gene)
    min_l = t.get("min_len", DEFAULTS["min_len"])
    max_l = t.get("max_len", DEFAULTS["max_len"])
    print(f"  [1/5] Descarregando sequências do NCBI (até {max_s})...")

    def _fetch_in_batches(ids: list[str], batch: int = 200) -> list[SeqRecord]:
        """efetch em lotes — necessário para centenas de IDs (evita URLs longos)."""
        recs: list[SeqRecord] = []
        for i in range(0, len(ids), batch):
            chunk = ids[i:i + batch]
            time.sleep(0.4)                     # respeitar limite NCBI (3 req/s sem api_key)
            h   = Entrez.efetch(db="nucleotide", id=",".join(chunk),
                                rettype="fasta", retmode="text")
            raw = h.read()
            h.close()
            recs.extend(_parse_efetch_fasta(raw))
        return recs

    def _run_query(query: str) -> list[SeqRecord]:
        for attempt in range(3):
            try:
                time.sleep(0.5)
                h   = Entrez.esearch(db="nucleotide", term=query, retmax=max_s)
                ids = Entrez.read(h)["IdList"]; h.close()
                if not ids:
                    return []
                recs = _fetch_in_batches(ids)
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
                raw = h.read()
                h.close()
                records = _parse_efetch_fasta(raw)
            except Exception as e:
                raise RuntimeError(f"Sem dados para {gene_key}: {e}")
        else:
            raise RuntimeError(f"Sem dados para {gene_key}")

    print(f"    ✔ {len(records)} sequências  |  modo: alinhamento múltiplo")
    return records

def select_length_cluster(records: list[SeqRecord], tol: float):
    """Mantém só o cluster de comprimento dominante (banda ±tol em torno do
    comprimento mais denso). Remove automaticamente fragmentos parciais e registos
    sobredimensionados que tornariam o alinhamento gappy — substitui o ajuste manual
    de min_len/max_len por gene, por isso funciona out-of-the-box em genes novos.

    Devolve (records_filtrados, (lo, hi, n_removidas)) ou (records, None) se não aplicado.
    """
    if tol <= 0 or len(records) < 4:
        return records, None
    lens = [len(r.seq) for r in records]
    # centro = comprimento com mais vizinhos dentro de ±tol (robusto a bimodalidade)
    best_c, best_n = lens[0], -1
    for c in set(lens):
        lo, hi = c * (1 - tol), c * (1 + tol)
        n = sum(lo <= L <= hi for L in lens)
        if n > best_n:
            best_c, best_n = c, n
    lo, hi = best_c * (1 - tol), best_c * (1 + tol)
    kept = [r for r in records if lo <= len(r.seq) <= hi]
    if len(kept) < 2:                       # salvaguarda: não reduzir abaixo do alinhável
        return records, None
    return kept, (int(lo), int(hi), len(records) - len(kept))

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
    fa_out.write_text(result.stdout, encoding="utf-8")
    aln = AlignIO.read(fa_out, "fasta")
    print(f"    ✔ {len(aln)} seqs × {aln.get_alignment_length()} posições")
    return aln

# ── 3. Janelas candidatas ────────────────────────────────────────────────────
# Conservação = PPI (Percentage of Pairwise Identity) — portado do ViruScope
# (Lima et al., anasfplima/ViruScope, viruscopeCLI.py). Em vez da frequência da
# base mais comum por coluna, usa-se a identidade par-a-par média (mais rigorosa),
# com pesos IUPAC para bases ambíguas. Decisão do utilizador: só PPI (sem PPI3).
_AMBIGUITY_MAP = {
    "A": {"A": 1.0, "R": 0.5, "W": 0.5, "M": 0.5, "D": 0.33, "H": 0.33, "V": 0.33, "N": 0.25},
    "C": {"C": 1.0, "Y": 0.5, "M": 0.5, "S": 0.5, "B": 0.33, "H": 0.33, "V": 0.33, "N": 0.25},
    "G": {"G": 1.0, "R": 0.5, "K": 0.5, "S": 0.5, "B": 0.33, "D": 0.33, "V": 0.33, "N": 0.25},
    "T": {"T": 1.0, "Y": 0.5, "K": 0.5, "W": 0.5, "B": 0.33, "D": 0.33, "H": 0.33, "N": 0.25},
    "R": {"A": 0.5, "G": 0.5, "R": 1.0, "N": 0.25},
    "Y": {"C": 0.5, "T": 0.5, "Y": 1.0, "N": 0.25},
    "S": {"C": 0.5, "G": 0.5, "S": 1.0, "N": 0.25},
    "W": {"A": 0.5, "T": 0.5, "W": 1.0, "N": 0.25},
    "K": {"G": 0.5, "T": 0.5, "K": 1.0, "N": 0.25},
    "M": {"A": 0.5, "C": 0.5, "M": 1.0, "N": 0.25},
    "B": {"C": 0.33, "G": 0.33, "T": 0.33, "B": 1.0, "N": 0.25},
    "D": {"A": 0.33, "G": 0.33, "T": 0.33, "D": 1.0, "N": 0.25},
    "H": {"A": 0.33, "C": 0.33, "T": 0.33, "H": 1.0, "N": 0.25},
    "V": {"A": 0.33, "C": 0.33, "G": 0.33, "V": 1.0, "N": 0.25},
    "N": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25, "N": 1.0},
}

def _pairwise_identity(col: list[str]) -> float:
    """Identidade par-a-par de uma coluna do alinhamento (ViruScope). Fração 0–1.

    Colunas com <30% de bases não-gap devolvem 0.0; pares de bases ambíguas
    contribuem com peso parcial via _AMBIGUITY_MAP.
    """
    column = [b.upper() for b in col]
    nogaps = [b for b in column if b != "-"]
    if len(nogaps) < 0.3 * len(column):
        return 0.0
    counts: dict[str, int] = {}
    for b in nogaps:
        counts[b] = counts.get(b, 0) + 1
    n = len(nogaps)
    total_pairs = n * (n - 1) / 2
    if total_pairs <= 0:
        return 0.0
    identical = sum(c * (c - 1) / 2 for c in counts.values())
    bases = list(counts)
    for i in range(len(bases)):
        for j in range(i + 1, len(bases)):
            b1, b2 = bases[i], bases[j]
            w = _AMBIGUITY_MAP.get(b1, {}).get(b2)
            if w:
                identical += counts[b1] * counts[b2] * w
    return identical / total_pairs

# Nota: a PPI por coluna é independente da janela, por isso é pré-calculada uma
# única vez por alinhamento em candidate_windows (essencial para centenas de seqs).

def candidate_windows(aln: AlignIO.MultipleSeqAlignment, gene_key: str) -> list[dict]:
    cons_min = cfg(gene_key, "cons_min")
    gap_max  = cfg(gene_key, "gap_max")
    plen_min = cfg(gene_key, "probe_len_min")
    plen_max = cfg(gene_key, "probe_len_max")
    step     = cfg(gene_key, "window_step")
    aln_len  = aln.get_alignment_length()
    n_seqs   = len(aln)
    cols     = [[str(rec.seq[i]).upper() for rec in aln] for i in range(aln_len)]

    # Pré-cálculo por coluna (uma vez): gap, validade, PPI e consenso
    col_gap   = [c.count("-") / n_seqs for c in cols]
    col_valid = [sum(1 for b in c if b != "-") >= 2 for c in cols]
    col_ppi   = [(_pairwise_identity(c) if v else 0.0) for c, v in zip(cols, col_valid)]
    col_cons  = [max(set(c) - {"-"}, key=c.count, default="N") for c in cols]

    windows = []
    for start in range(0, aln_len - plen_min + 1, step):
        for plen in range(plen_min, plen_max + 1):
            end = start + plen
            if end > aln_len:
                break
            if any(col_gap[i] > gap_max for i in range(start, end)):
                continue
            if n_seqs < 2:                       # seq única → PPI indefinida, cons=1.0
                cons = 1.0
            else:
                vals = [col_ppi[i] for i in range(start, end) if col_valid[i]]
                cons = (sum(vals) / len(vals)) if vals else 0.0
            if cons < cons_min:
                continue
            consensus = "".join(col_cons[i] for i in range(start, end))
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

# No-fold score (ViruScope): converte ΔG (kcal/mol) num score 0–100 via penalidade
# logística — alternativa contínua ao pass/fail rígido. Score alto = pouca estrutura
# secundária indesejada. Mantemos também os flags pass/fail para compatibilidade.
def fold_score(delta_g_kcal: float, dg_threshold_cal: float = -7000.0,
               slope: float = 1.0, T: float = 310.0) -> Optional[float]:
    if delta_g_kcal is None or math.isnan(delta_g_kcal):
        return None
    R = 1.987  # cal/(mol·K)
    delta_g_cal = delta_g_kcal * 1000.0  # ViruScope trabalha em cal/mol
    x = (delta_g_cal - (dg_threshold_cal / 2)) / (slope * R * T)
    x = max(-700.0, min(700.0, x))       # clamp → evita OverflowError (ΔG enorme/±inf)
    penalty = 1.0 / (1.0 + math.exp(x))
    return round((1.0 - penalty) * 100.0, 2)

def nofold_score(probe: "Probe") -> Optional[float]:
    """Média dos fold_scores de hairpin, homodímero e seqfold ΔG (os disponíveis)."""
    scores = [fold_score(dg) for dg in (probe.hairpin_dg, probe.homodimer_dg,
                                        probe.seqfold_dg) if dg is not None]
    scores = [s for s in scores if s is not None]
    return round(sum(scores) / len(scores), 2) if scores else None

def probe_quality(probe: "Probe") -> float:
    """Score de qualidade 0–1 (conservação PPI + No-fold, peso igual) para priorizar
    quais probes vão ao Boltz — gasta o tempo de GPU nas melhores, não nas de Tm mais alto.
    Tm/GC/estrutura já foram filtrados pelos passes; isto ordena os sobreviventes."""
    cons = probe.conservation if probe.conservation is not None else 0.0
    nof  = (probe.nofold_score / 100.0) if probe.nofold_score is not None else 0.0
    return (cons + nof) / 2

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
        dg = sf.dg(seq, temp=37.0)
        if not math.isfinite(dg):        # seqfold pode devolver +-inf p/ certas seqs
            probe.notes += "[seqfold dg nao-finito (ignorado)] "
            return probe                  # seqfold_dg fica None
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
    for _old in yaml_dir.glob("*.yaml"):     # limpar YAMLs de corridas anteriores
        _old.unlink()                         # (evita acumular probes desatualizadas no zip)

    selected: list[Probe] = []
    for gk in TARGETS:
        cands = sorted(
            [p for p in probes if p.gene == gk and p.pass_basic and p.pass_seqfold],
            key=probe_quality, reverse=True            # melhores por qualidade (PPI + No-fold)
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

    meta_fields = ["probe_id", "gene", "organism", "sequence", "tm", "gc",
                   "conservation", "hairpin_dg", "homodimer_dg", "seqfold_dg", "nofold_score"]
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
    print(f"  {len(selected)} probes seleccionadas (top {top_n}/gene por qualidade: PPI + No-fold)")
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

def export_colab_from_csv(top_n: int) -> Path:
    """Gera os inputs do Colab a partir do FINAL_PROBES_ALL.csv já existente
    (sem re-correr NCBI/MAFFT) — determinístico e consistente com o CSV commitado.
    Reutiliza export_colab_inputs(). Produz output/colab_boltz2/boltz2_inputs.zip
    pronto para upload direto no notebook colab_boltz2_batch.ipynb."""
    csv_path = BASE_DIR / "output" / "FINAL_PROBES_ALL.csv"
    df = pd.read_csv(csv_path)
    df = df[df["fonte"] == "romeu"]

    def _f(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    def _b(v):
        return str(v).strip().lower() in ("true", "1", "1.0")

    probes = [
        Probe(
            gene=r["gene"], organism=r["organism"], group=str(r.get("group", "")),
            probe_id=r["probe_id"], sequence=str(r["sequence"]),
            tm=_f(r["tm"]) or 0.0, gc=_f(r["gc"]) or 0.0,
            conservation=_f(r.get("conservation")), nofold_score=_f(r.get("nofold_score")),
            hairpin_dg=_f(r["hairpin_dg"]), homodimer_dg=_f(r["homodimer_dg"]),
            seqfold_dg=_f(r["seqfold_dg"]),
            pass_basic=_b(r["pass_basic"]), pass_seqfold=_b(r["pass_seqfold"]),
        )
        for _, r in df.iterrows()
    ]
    print(f"  {len(probes)} probes próprias lidas de {csv_path.name}")
    return export_colab_inputs(probes, top_n=top_n)

# ── 7. Outputs por gene ──────────────────────────────────────────────────────
def write_gene_outputs(probes: list[Probe], gene_key: str):
    out_dir = ALIGN_DIR / gene_key
    tsv   = out_dir / f"{gene_key}_probes_scored.tsv"
    fasta = out_dir / f"{gene_key}_viroscope_probes.fasta"
    fields = list(asdict(Probe("", "", "")).keys())   # robusto a lista vazia
    with open(tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for p in probes:
            w.writerow(asdict(p))
    with open(fasta, "w", encoding="utf-8") as f:
        for p in probes:
            if p.pass_basic:
                f.write(f">{p.probe_id}\n{p.sequence}\n")
    print(f"    ✔ {tsv.name}  ✔ {fasta.name}  ({len(probes)} probes)")

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

    try:
        xl = pd.ExcelFile(path)
    except ImportError as e:
        print(f"  ⚠ openpyxl em falta — probes da Beatriz ignoradas "
              f"(instala com: pip install openpyxl). [{e}]")
        return rows

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
    # Janela Tm derivada da temperatura do ensaio, se definida
    if DEFAULTS.get("assay_temp_c") is not None:
        at = DEFAULTS["assay_temp_c"]
        DEFAULTS["tm_min"] = round(at + DEFAULTS["tm_offset_min"], 1)
        DEFAULTS["tm_max"] = round(at + DEFAULTS["tm_offset_max"], 1)
        print(f"  Ensaio @ {at}°C → janela Tm derivada: {DEFAULTS['tm_min']}–{DEFAULTS['tm_max']}°C "
              f"(overrides AT-rich por gene mantêm-se)")
    print("═"*60)

    all_probes: list[Probe] = []
    probe_counter = itertools.count(1)

    for gene_key, t in TARGETS.items():
        print(f"\n{'━'*60}")
        print(f"  {gene_key.upper()}  |  {t['organism']}  |  Grupo {t['group']}")
        print(f"{'━'*60}")

        records = fetch_sequences(gene_key)
        if cfg(gene_key, "len_cluster"):
            records, info = select_length_cluster(records, cfg(gene_key, "len_cluster_tol"))
            if info:
                lo, hi, dropped = info
                print(f"  [1b] Cluster de comprimento dominante: {lo}–{hi} bp "
                      f"→ {len(records)} seqs (removidas {dropped} fora do cluster)")
        aln     = align_mafft(records, gene_key)
        windows = candidate_windows(aln, gene_key)

        probes: list[Probe] = []
        print(f"  [4/5] Scoring básico (Tm, GC, hairpin, homodimer)...")
        for w in windows:
            sc = score_probe(w["seq"])
            p  = Probe(
                gene=gene_key, organism=t["organism"], group=t["group"],
                sequence=w["seq"], pos_start=w["start"], pos_end=w["end"],
                tm=sc["tm"], gc=sc["gc"], conservation=round(w["cons"], 3),
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

        # No-fold score (ViruScope) — depende de seqfold_dg, por isso após o seqfold
        for p in probes:
            p.nofold_score = nofold_score(p)

        write_gene_outputs(probes, gene_key)
        all_probes.extend(probes)

    write_consolidated_csv(all_probes)

    # ── Resumo por gene ───────────────────────────────────────────────────────
    print("\n" + "═"*60)
    print("  RESUMO POR GENE")
    print("═"*60)
    print(f"  {'Gene':<8} {'Janelas':>8} {'Pass Tm':>8} {'Pass GC':>8} {'Pass HP':>8} {'Pass Dim':>9} {'Pass Básico':>11} {'Pass Seqfold':>12}")
    print(f"  {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*9} {'─'*11} {'─'*12}")

    totals = [0] * 7
    for gk in TARGETS:
        ps = [p for p in all_probes if p.gene == gk]
        if not ps:
            continue
        n_tot = len(ps)
        tm_min = cfg(gk,"tm_min"); tm_max = cfg(gk,"tm_max")
        gc_min = cfg(gk,"gc_min"); gc_max = cfg(gk,"gc_max")
        hp_min = cfg(gk,"hp_min"); dm_min = cfg(gk,"dimer_min")
        n_tm  = sum(1 for p in ps if tm_min <= p.tm  <= tm_max)
        n_gc  = sum(1 for p in ps if gc_min <= p.gc  <= gc_max and tm_min <= p.tm <= tm_max)
        n_hp  = sum(1 for p in ps if gc_min <= p.gc  <= gc_max and tm_min <= p.tm <= tm_max
                    and (p.hairpin_dg is None or p.hairpin_dg >= hp_min))
        n_bas = sum(p.pass_basic   for p in ps)
        n_sf  = sum(p.pass_seqfold for p in ps)
        for i, v in enumerate([n_tot, n_tm, n_gc, n_hp, n_bas, n_bas, n_sf]):
            totals[i] += v
        print(f"  {gk:<8} {n_tot:>8} {n_tm:>8} {n_gc:>8} {n_hp:>8} {n_bas:>9} {n_bas:>11} {n_sf:>12}")

    print(f"  {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*9} {'─'*11} {'─'*12}")
    print(f"  {'TOTAL':<8} {totals[0]:>8} {totals[1]:>8} {totals[2]:>8} {totals[3]:>8} {totals[4]:>9} {totals[4]:>11} {totals[6]:>12}")

    # ── Funil global ──────────────────────────────────────────────────────────
    print(f"\n  FUNIL GLOBAL (todos os genes)")
    print(f"  {'─'*50}")
    total = len(all_probes)
    remaining = all_probes[:]
    steps_global = [
        ("Janelas candidatas totais", None),
        ("Pass Tm",       lambda p: cfg(p.gene,"tm_min")  <= p.tm <= cfg(p.gene,"tm_max")),
        ("Pass GC",       lambda p: cfg(p.gene,"gc_min")  <= p.gc <= cfg(p.gene,"gc_max")),
        ("Pass Hairpin",  lambda p: p.hairpin_dg  is None or p.hairpin_dg  >= cfg(p.gene,"hp_min")),
        ("Pass Homodimer",lambda p: p.homodimer_dg is None or p.homodimer_dg >= cfg(p.gene,"dimer_min")),
    ]
    current = total
    print(f"  {'Janelas candidatas totais':<30} {total:>6}")
    for label, test in steps_global[1:]:
        ok   = [p for p in remaining if test(p)]
        fail = len(remaining) - len(ok)
        print(f"  ↳ {label:<28} {len(ok):>6}  [-{fail}]")
        remaining = ok
    n_sf_total = sum(p.pass_seqfold for p in all_probes)
    print(f"  ↳ {'Pass seqfold':<28} {n_sf_total:>6}  [-{len(remaining)-n_sf_total}]")
    print(f"  {'─'*50}")
    print(f"  {'Aprovadas para Colab/Boltz-2':<30} {n_sf_total:>6}")

    if colab_top > 0:
        export_colab_inputs(all_probes, top_n=colab_top)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="GFET Probe Pipeline")
    ap.add_argument("--no-seqfold", action="store_true",
                    help="Saltar análise de estrutura secundária (seqfold)")
    ap.add_argument("--colab", type=int, default=0, metavar="N",
                    help="Gerar YAMLs + ZIP para Boltz-2 Colab (top N por gene, ex: --colab 5)")
    ap.add_argument("--max-seqs", type=int, default=None, metavar="N",
                    help=f"Sequências NCBI a descarregar por gene (default {DEFAULTS['max_seqs']}). "
                         f"Escalar com rigor — ver scripts/analysis.py (rarefação).")
    ap.add_argument("--no-cluster", action="store_true",
                    help="Desativar a seleção automática do cluster de comprimento dominante")
    ap.add_argument("--cluster-tol", type=float, default=None, metavar="T",
                    help=f"Tolerância do cluster de comprimento, ±T (default {DEFAULTS['len_cluster_tol']})")
    ap.add_argument("--assay-temp", type=float, default=None, metavar="C",
                    help="Temperatura do ensaio de hibridação (°C) → deriva a janela Tm "
                         "automaticamente (T+15 a T+35). Omitir = manter janela default 53–72.")
    ap.add_argument("--export-colab", type=int, default=0, metavar="N",
                    help="Gerar inputs Colab (top N/gene) a partir do FINAL_PROBES_ALL.csv "
                         "existente, SEM re-correr o pipeline. Produz boltz2_inputs.zip.")
    args = ap.parse_args()
    if args.export_colab > 0:                       # atalho: só (re)gerar inputs Colab
        export_colab_from_csv(args.export_colab)
        raise SystemExit(0)
    if args.max_seqs is not None:
        DEFAULTS["max_seqs"] = args.max_seqs
    if args.assay_temp is not None:
        DEFAULTS["assay_temp_c"] = args.assay_temp
    if args.no_cluster:
        DEFAULTS["len_cluster"] = False
    if args.cluster_tol is not None:
        DEFAULTS["len_cluster_tol"] = args.cluster_tol
    run_pipeline(run_seqfold=not args.no_seqfold, colab_top=args.colab)
