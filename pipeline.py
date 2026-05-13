"""
GFET Probe Pipeline
─────────────────────────────────────────────────────────────────────────────
Pipeline computacional para design in silico de probes DNA para biossensores
GFET (Graphene Field-Effect Transistor) em detecção de patogénios bacterianos.

Etapas:
  NCBI → MAFFT → janelas conservadas → primer3 → seqfold → Boltz-2 (3D)

Referências dos limiares:
  SantaLucia & Hicks 2004  — Tm nearest-neighbour, parâmetros primer3
  Wetmur 1991              — comprimento da probe (18–28 nt)
  Zadeh et al. 2011        — limiar seqfold (proxy NUPACK ensemble)
  IDT OligoAnalyzer        — GC 40–65%, homodimer ΔG > -5 kcal/mol
  Stover et al. 2000       — P. aeruginosa GC genómico ~67% (PAO1)
  Wohlwend et al. 2024     — Boltz-2, predição estrutura 3D
─────────────────────────────────────────────────────────────────────────────
"""

import os, sys, re, csv, json, yaml, shutil, subprocess, tempfile, textwrap, itertools, time
from io import StringIO
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional

import requests
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO
from Bio.PDB import PDBParser, Superimposer, MMCIFIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import primer3

# ─── Entrez ────────────────────────────────────────────────────────────────
Entrez.email = "pg45861@uminho.pt"

# ─── Thermodynamic parameters (primer3) ────────────────────────────────────
# SantaLucia & Hicks 2004 — nearest-neighbour, condições fisiológicas/GFET
P3_MV_CONC  = 50.0     # [Na+] mM
P3_DV_CONC  = 0.0      # [Mg2+] mM
P3_DNTP     = 0.0      # [dNTP] mM
P3_DNA_CONC = 250.0    # [oligo] nM
P3_TEMP_C   = 37.0     # temperatura de referência

# ─── Paths ─────────────────────────────────────────────────────────────────
BASE_DIR     = Path(__file__).parent
ALIGN_DIR    = BASE_DIR / "output" / "alignments"
STRUCT_DIR   = BASE_DIR / "output" / "structures"
BEATRIZ_XLSX = BASE_DIR / "data" / "BeatrizMasterThesis_ProbesData_IPLEXMED_ENERGIAS_CAPS.xlsx"

# MAFFT: procura no projeto primeiro (duas estruturas possíveis), depois PATH, depois Tese
_MAFFT_CANDIDATES = [
    str(BASE_DIR / "MAFFT" / "mafft-win" / "mafft.bat"),
    str(BASE_DIR / "MAFFT" / "mafft-7.526-win64-signed" / "mafft-win" / "mafft.bat"),
    shutil.which("mafft"),
    shutil.which("mafft.bat"),
    str(Path.home() / "OneDrive" / "Ambiente de Trabalho" / "Tese"
        / "MAFFT" / "mafft-7.526-win64-signed" / "mafft-win" / "mafft.bat"),
]
MAFFT_BIN = Path(next(p for p in _MAFFT_CANDIDATES if p and Path(p).exists()))

# NucleoFold3D: actualizar se o script estiver noutro sítio
NUCLEOFOLD3D_SCRIPT = (
    Path.home() / "OneDrive" / "Ambiente de Trabalho" / "Tese"
    / "NucleoFold" / "NucleoFold3D.py"
)

for d in [ALIGN_DIR, STRUCT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ─── Targets ───────────────────────────────────────────────────────────────
# Cada entrada pode sobrepor limiares globais com chaves opcionais:
#   tm_min, tm_max, gc_min, gc_max, hp_min, dimer_min, cons_min, gap_max,
#   probe_len_min, probe_len_max, nupack_dg_max, nupack_defect_max
TARGETS = {
    "nuc": {
        "organism":          "Staphylococcus aureus",
        "group":             "A",
        "ncbi_query":        'nuc[Gene Name] AND "Staphylococcus aureus"[Organism] AND 200:3000[Sequence Length]',
        "ncbi_query_alt":    '"Staphylococcus aureus"[Organism] AND nuc[Title] AND 200:2000[Sequence Length]',
        "refseq_fallback":   "NC_002951",
        "min_len":           200,
        "max_len":           3000,
        "tm_min":            52.0,  # gene AT-rico — SantaLucia & Hicks 2004
        "gc_min":            0.38,  # AT-rico, abaixo do padrão 40% — IDT OligoAnalyzer
        "gc_max":            0.60,  # AT-rico, limitar máximo — IDT OligoAnalyzer
    },
    "rmpM": {
        "organism":          "Neisseria meningitidis",
        "group":             "A",
        "ncbi_query":        'rmpM[Gene Name] AND "Neisseria meningitidis"[Organism] AND 100:1200[Sequence Length]',
        "ncbi_query_alt":    '("Neisseria meningitidis"[Organism]) AND (rmpM[Title] OR "class 4 outer membrane protein"[Title]) AND 100:1200[Sequence Length]',
        "ncbi_query_fallback": '"Neisseria meningitidis"[Organism] AND "outer membrane protein"[Title] AND 100:1200[Sequence Length]',
        "refseq_fallback":   "NC_003112",
        "min_len":           100,
        "max_len":           1200,
        "gc_max":            0.65,  # IDT OligoAnalyzer guidelines
        "cons_min":          0.80,  # limiar permissivo — poucos homólogos NCBI
        "allow_single_seq":  True,  # fallback para sequência única se <2 homólogos
    },
    "lytA": {
        "organism":          "Streptococcus pneumoniae",
        "group":             "B",
        "ncbi_query":        'lytA[Gene Name] AND "Streptococcus pneumoniae"[Organism] AND 700:1300[Sequence Length]',
        "ncbi_query_alt":    '"Streptococcus pneumoniae"[Organism] AND lytA[Title] AND 700:1300[Sequence Length]',
        "refseq_fallback":   "NC_003028",
        "min_len":           700,
        "max_len":           1300,
        "cons_min":          0.70,  # lytA tem diversidade alélica — Whatmore et al. 2000
        "gap_max":           0.40,
    },
    "oprL": {
        "organism":          "Pseudomonas aeruginosa",
        "group":             "B",
        "ncbi_query":        'oprL[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 300:2000[Sequence Length]',
        "ncbi_query_alt":    '"Pseudomonas aeruginosa"[Organism] AND oprL[Title] AND 300:2000[Sequence Length]',
        "refseq_fallback":   "NC_002516",
        "min_len":           300,
        "max_len":           2000,
        "gc_max":            0.70,  # P. aeruginosa ~67% GC genómico — Stover et al. 2000
    },
    "algD": {
        "organism":          "Pseudomonas aeruginosa",
        "group":             "B",
        "ncbi_query":        'algD[Gene Name] AND "Pseudomonas aeruginosa"[Organism] AND 500:2500[Sequence Length]',
        "ncbi_query_alt":    '"Pseudomonas aeruginosa"[Organism] AND algD[Title] AND 500:2500[Sequence Length]',
        "refseq_fallback":   "NC_002516",
        "min_len":           500,
        "max_len":           2500,
        "gc_max":            0.70,  # P. aeruginosa ~67% GC genómico — Stover et al. 2000
    },
    "frdB": {
        "organism":          "Haemophilus influenzae",
        "group":             "B",
        "ncbi_query":        'frdB[Gene Name] AND "Haemophilus influenzae"[Organism] AND 200:2000[Sequence Length]',
        "ncbi_query_alt":    '"Haemophilus influenzae"[Organism] AND frdB[Title] AND 200:2000[Sequence Length]',
        "refseq_fallback":   "NC_000907",
        "min_len":           200,
        "max_len":           2000,
        "gc_min":            0.38,  # H. influenzae AT-rico — IDT OligoAnalyzer
        "gc_max":            0.60,
    },
}

# ─── Limiares globais (podem ser sobrepostos por gene em TARGETS) ───────────
DEFAULTS = {
    "tm_min":             53.0,   # °C  — SantaLucia & Hicks 2004
    "tm_max":             72.0,   # °C
    "gc_min":             0.40,   # IDT OligoAnalyzer guidelines
    "gc_max":             0.65,
    "hp_min":            -2.0,    # kcal/mol — hairpin ΔG mínimo (ssDNA, Wetmur 1991)
    "dimer_min":         -5.0,    # kcal/mol — homodimer ΔG mínimo (IDT OligoAnalyzer)
    "cons_min":           0.85,   # conservação mínima — Zadeh et al. 2011
    "gap_max":            0.20,   # gap máximo no alinhamento
    "probe_len_min":      18,     # nt — Wetmur 1991
    "probe_len_max":      28,     # nt
    "window_step":         3,     # passo da janela deslizante
    "max_seqs":           20,     # seqs a descarregar do NCBI
    "min_len":           100,     # bp — comprimento mínimo da sequência alvo
    "max_len":          5000,     # bp
    "nupack_dg_max":     -6.0,    # kcal/mol — ΔG MFE seqfold (proxy NUPACK, Zadeh 2011)
                                  # Rejeita probes com estrutura secundária forte (< -6 kcal/mol)
    "nupack_defect_max":  0.50,   # fracção de bases emparelhadas no MFE (informativo)
    "rmsd_max":           2.0,    # Å — threshold de qualidade estrutural
}

def cfg(gene_key: str, param: str):
    return TARGETS[gene_key].get(param, DEFAULTS[param])

# ─── Dataclass probe ────────────────────────────────────────────────────────
@dataclass
class Probe:
    gene:          str
    organism:      str
    group:         str
    probe_id:      str = ""
    sequence:      str = ""
    pos_start:     int = 0
    pos_end:       int = 0
    tm:            float = 0.0
    gc:            float = 0.0
    hairpin_dg:    Optional[float] = None
    homodimer_dg:  Optional[float] = None
    nupack_dg:     Optional[float] = None
    nupack_defect: Optional[float] = None
    structure_cif: Optional[str] = None
    rmsd:          Optional[float] = None
    pass_basic:    bool = False
    pass_nupack:   bool = False
    pass_3d:       bool = False
    fonte:         str = "romeu"
    notes:         str = ""

# ─── 1. Nomenclatura ────────────────────────────────────────────────────────
_ORG_ABBR = {
    "Staphylococcus aureus":      "Saur",
    "Neisseria meningitidis":     "Nmen",
    "Streptococcus pneumoniae":   "Spne",
    "Pseudomonas aeruginosa":     "Paer",
    "Haemophilus influenzae":     "Hinf",
}

def build_probe_id(probe: Probe, index: int) -> str:
    """
    p001_Saur_nuc_pos373-395_Tm54.1_GC40_hp-0.38_PASS
    Formato: p{idx}_{org}_{gene}_pos{start}-{end}_Tm{tm}_GC{gc}_hp{hp}_{status}
    """
    org   = _ORG_ABBR.get(probe.organism, probe.organism[:4].capitalize())
    tm    = f"{probe.tm:.1f}"
    gc    = f"{int(round(probe.gc * 100))}"
    hp    = f"{probe.hairpin_dg:.2f}" if probe.hairpin_dg is not None else "NA"
    stat  = "PASS" if probe.pass_basic else "FAIL"
    return f"p{index:03d}_{org}_{probe.gene}_pos{probe.pos_start}-{probe.pos_end}_Tm{tm}_GC{gc}_hp{hp}_{stat}"

# ─── 2. NCBI download ────────────────────────────────────────────────────────
def fetch_sequences(gene_key: str) -> list[SeqRecord]:
    t     = TARGETS[gene_key]
    max_s = DEFAULTS["max_seqs"]
    min_l = t.get("min_len", DEFAULTS["min_len"])
    max_l = t.get("max_len", DEFAULTS["max_len"])
    print(f"  [1/6] Descarregando sequências do NCBI...")

    def _run_query(query: str) -> list[SeqRecord]:
        for attempt in range(3):
            try:
                time.sleep(0.5)
                h    = Entrez.esearch(db="nucleotide", term=query, retmax=max_s)
                ids  = Entrez.read(h)["IdList"]; h.close()
                if not ids:
                    return []
                time.sleep(0.5)
                h    = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
                raw  = "\n".join(l for l in h.read().splitlines()
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

    # allow_single_seq: genes com poucos homólogos (ex: rmpM)
    if len(records) < 2 and t.get("allow_single_seq") and len(records) == 1:
        print(f"    ⚠ Apenas 1 seq — modo sequência única ({gene_key}).")
        print(f"    ✔ 1 sequência  |  modo: sequência única")
        return records

    if not records:
        fb = t.get("refseq_fallback")
        if fb:
            print(f"    ⚠ Sem resultados — a usar referência {fb}")
            try:
                time.sleep(0.5)
                h    = Entrez.efetch(db="nucleotide", id=fb, rettype="fasta", retmode="text")
                raw  = "\n".join(l for l in h.read().splitlines()
                                 if not l.startswith(("!", "#", ";")))
                h.close()
                records = list(SeqIO.parse(StringIO(raw), "fasta"))
                print(f"    ✔ Referência {fb} descarregada.")
            except Exception as e:
                raise RuntimeError(f"Sem dados para {gene_key}: {e}")
        else:
            raise RuntimeError(f"Sem dados para {gene_key}")

    print(f"    ✔ {len(records)} sequências  |  modo: alinhamento múltiplo")
    return records

# ─── 3. MAFFT ────────────────────────────────────────────────────────────────
def align_mafft(records: list[SeqRecord], gene_key: str) -> AlignIO.MultipleSeqAlignment:
    out_dir = ALIGN_DIR / gene_key
    out_dir.mkdir(exist_ok=True)
    fa_in  = out_dir / "input.fasta"
    fa_out = out_dir / "aligned.fasta"
    if len(records) == 1:
        print(f"  [2/6] Alinhamento ignorado (sequência única).")
        SeqIO.write(records, fa_out, "fasta")
        return AlignIO.read(fa_out, "fasta")
    SeqIO.write(records, fa_in, "fasta")
    print(f"  [2/6] Alinhamento MAFFT ({len(records)} sequências)...")
    cmd = [str(MAFFT_BIN), "--auto", str(fa_in)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    fa_out.write_text(result.stdout)
    aln = AlignIO.read(fa_out, "fasta")
    print(f"    ✔ {len(aln)} seqs × {aln.get_alignment_length()} posições")
    return aln

# ─── 4. Janelas candidatas ────────────────────────────────────────────────────
def _conservation(col: list[str]) -> float:
    bases = [b for b in col if b not in ("-", "N", "n")]
    if not bases:
        return 0.0
    return max(bases.count(b) for b in set(bases)) / len(bases)

def _gap_freq(col: list[str]) -> float:
    return col.count("-") / len(col)

def candidate_windows(aln: AlignIO.MultipleSeqAlignment, gene_key: str) -> list[dict]:
    cons_min  = cfg(gene_key, "cons_min")
    gap_max   = cfg(gene_key, "gap_max")
    plen_min  = cfg(gene_key, "probe_len_min")
    plen_max  = cfg(gene_key, "probe_len_max")
    step      = cfg(gene_key, "window_step")
    aln_len   = aln.get_alignment_length()
    cols      = [[str(rec.seq[i]) for rec in aln] for i in range(aln_len)]
    windows   = []
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
            consensus = "".join(max(set(c) - {"-"}, key=lambda b: c.count(b), default="N")
                                 for c in w_cols).upper()
            windows.append({"start": start, "end": end, "cons": cons, "seq": consensus})
    print(f"  [3/6] ✔ {len(windows)} janelas candidatas (cons ≥ {cons_min}, gap ≤ {gap_max*100:.0f}%)")
    return windows

# ─── 5. Scoring básico (Tm, GC, hairpin, homodimer) ──────────────────────────
def score_probe(seq: str) -> dict:
    kw  = dict(mv_conc=P3_MV_CONC, dv_conc=P3_DV_CONC,
               dntp_conc=P3_DNTP, dna_conc=P3_DNA_CONC)
    tm  = primer3.calc_tm(seq, **kw)
    gc  = (seq.count("G") + seq.count("C")) / len(seq)
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
    ok = (
        cfg(gene_key, "tm_min")    <= probe.tm <= cfg(gene_key, "tm_max") and
        cfg(gene_key, "gc_min")    <= probe.gc <= cfg(gene_key, "gc_max") and
        (probe.hairpin_dg  is None or probe.hairpin_dg  >= cfg(gene_key, "hp_min")) and
        (probe.homodimer_dg is None or probe.homodimer_dg >= cfg(gene_key, "dimer_min"))
    )
    return ok

# ─── 6. seqfold (substituto open-source do NUPACK) ───────────────────────────
def run_seqfold_probe(probe: Probe, gene_key: str) -> Probe:
    """
    Análise de estrutura secundária via seqfold (Nussinov, MIT license).
    Critério de aprovação: ΔG MFE ≥ nupack_dg_max (padrão -6.0 kcal/mol).
    nupack_defect = fracção de bases emparelhadas (calculado mas apenas informativo).
    Ref: Zadeh et al. 2011 (NUPACK, referência conceptual); seqfold (pip install seqfold)
    """
    try:
        import seqfold as sf
    except ImportError:
        probe.notes += "[seqfold não instalado: pip install seqfold] "
        return probe

    seq = probe.sequence
    try:
        dg      = sf.dg(seq, temp=37.0)
        structs = sf.fold(seq, temp=37.0)
        # s.ij é lista de tuplas (i, j) — cada par representa uma ligação de bases
        paired  = {idx for s in structs for pair in s.ij for idx in pair if idx >= 0}
        paired_frac = len(paired) / len(seq) if seq else 0.0
        probe.nupack_dg     = round(dg, 2)
        probe.nupack_defect = round(paired_frac, 3)
    except Exception as e:
        probe.notes += f"[seqfold erro: {e}] "
        return probe

    probe.pass_nupack = probe.nupack_dg >= cfg(gene_key, "nupack_dg_max")
    return probe

# ─── 7. Estrutura 3D ─────────────────────────────────────────────────────────

def run_nucleofold(probe: Probe) -> Optional[Path]:
    """
    Corre o pipeline local NucleoFold3D.py para uma probe.
    Requer: nsp binary, RNAfold (ViennaRNA), AmberTools (tleap/sander/ambpdb).
    RNAfold e AmberTools funcionam em Windows via conda.
    O binário nsp é normalmente Linux — usar WSL se necessário.

    Devolve o path do PDB de menor energia ou None se falhar/não instalado.
    Ref: NucleoFold3D pipeline (Beatriz et al.) — RNAfold → NSP → AMBER OL15
    """
    if not NUCLEOFOLD3D_SCRIPT.exists():
        probe.notes += "[NucleoFold3D.py não encontrado] "
        return None

    # Criar CSV de input no formato esperado pelo NucleoFold3D.py
    tmp_csv = STRUCT_DIR / f"{probe.probe_id}_nf3d_input.csv"
    tmp_csv.write_text(f"job_name,sequence\n{probe.probe_id},{probe.sequence}\n",
                       encoding="utf-8")
    try:
        subprocess.run(
            [sys.executable, str(NUCLEOFOLD3D_SCRIPT), "--csv", str(tmp_csv)],
            check=True, capture_output=True, timeout=3600
        )
    except FileNotFoundError:
        probe.notes += "[NucleoFold3D: python não encontrado] "
        return None
    except subprocess.CalledProcessError as e:
        probe.notes += f"[NucleoFold3D erro: {e.returncode}] "
        return None
    except subprocess.TimeoutExpired:
        probe.notes += "[NucleoFold3D timeout] "
        return None
    finally:
        tmp_csv.unlink(missing_ok=True)

    # NucleoFold3D guarda em ~/scripts/results_<input_stem>/<job_name>/
    results_base = Path.home() / "scripts" / f"results_{tmp_csv.stem}"
    min_pdb = next(results_base.glob(f"{probe.probe_id}*/*_minenerg.pdb"), None)
    if min_pdb is None:
        # fallback: qualquer minenerg
        min_pdb = next(results_base.glob("**/*_minenerg.pdb"), None)
    if min_pdb:
        dest = STRUCT_DIR / f"{probe.probe_id}_nucleofold_minenerg.pdb"
        shutil.copy(min_pdb, dest)
        return dest
    return None

def _boltz_cmd() -> list[str] | None:
    """
    Devolve o comando para invocar o Boltz-2.
    Tenta 'boltz' no PATH primeiro; se não existir, usa 'python -m boltz'
    (instalação via pip sem entry-point no PATH).
    """
    if shutil.which("boltz"):
        return ["boltz"]
    try:
        import importlib.util
        if importlib.util.find_spec("boltz") is not None:
            return [sys.executable, "-m", "boltz"]
    except Exception:
        pass
    return None

def run_boltz2(probe: Probe, n_samples: int = 3) -> list[Path]:
    """
    Corre o Boltz-2 localmente para gerar n_samples estruturas CIF.
    Usa o CLI: boltz predict <yaml> --out_dir <dir> --diffusion_samples n
    Ref: Wohlwend et al. 2024 (Boltz-2 paper)
    """
    cmd_base = _boltz_cmd()
    if cmd_base is None:
        probe.notes += "[Boltz-2 não encontrado: pip install boltz] "
        return []
    out_dir = STRUCT_DIR / probe.probe_id / "boltz"
    out_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = out_dir / "input.yaml"
    yaml_path.write_text(yaml.dump({
        "version": 1,
        "sequences": [{"dna": {"id": ["A"], "sequence": probe.sequence}}]
    }))
    try:
        subprocess.run(
            cmd_base + [
                "predict", str(yaml_path),
                "--out_dir", str(out_dir),
                "--recycling_steps", "3",
                "--sampling_steps", "200",
                "--diffusion_samples", str(n_samples),
                "--accelerator", "cpu",
                "--model", "boltz2",
            ],
            check=True, capture_output=True, timeout=1800
        )
    except Exception as e:
        probe.notes += f"[Boltz erro: {e}] "
        return []
    return sorted(out_dir.glob("**/*.cif"))

# ─── 8. RMSD entre réplicas ───────────────────────────────────────────────────
def _calc_rmsd_between_cifs(cif_paths: list[Path]) -> Optional[float]:
    """
    Calcula RMSD médio entre pares de estruturas CIF usando Biopython Superimposer.
    Ref: Hamelryck & Manderick 2003 (Bio.PDB)
    Threshold de aprovação: RMSD_MAX (padrão 2.0 Å)
    """
    if len(cif_paths) < 2:
        return None
    parser    = PDBParser(QUIET=True)
    sup       = Superimposer()
    rmsds     = []
    for a, b in itertools.combinations(cif_paths, 2):
        try:
            s1 = parser.get_structure("s1", str(a))
            s2 = parser.get_structure("s2", str(b))
            atoms1 = [at for at in s1.get_atoms() if at.get_name() == "C4'"]
            atoms2 = [at for at in s2.get_atoms() if at.get_name() == "C4'"]
            n = min(len(atoms1), len(atoms2))
            if n < 3:
                continue
            sup.set_atoms(atoms1[:n], atoms2[:n])
            rmsds.append(sup.rms)
        except Exception:
            continue
    return sum(rmsds) / len(rmsds) if rmsds else None

def run_3d_pipeline(probe: Probe) -> Probe:
    """
    Tenta NucleoFold3D.py (local) primeiro; se falhar, tenta Boltz-2.
    NucleoFold3D gera 5 PDB minimizados — usa o de menor energia AMBER.
    Boltz-2 gera n réplicas CIF — calcula RMSD entre elas.
    Pass_3d = True se NucleoFold3D produziu resultado OU RMSD Boltz < rmsd_max.
    """
    print(f"      3D: NucleoFold3D... ", end="", flush=True)
    pdb = run_nucleofold(probe)
    if pdb:
        probe.structure_cif = str(pdb)   # campo reutilizado para PDB também
        probe.pass_3d       = True
        print(f"✔ ({pdb.name})")
        return probe

    print("falhou → Boltz-2... ", end="", flush=True)
    cifs = run_boltz2(probe, n_samples=3)
    if not cifs:
        print("✘")
        return probe
    probe.structure_cif = str(cifs[0])
    rmsd = _calc_rmsd_between_cifs(cifs)
    probe.rmsd  = rmsd
    rmsd_thresh = DEFAULTS["rmsd_max"]
    probe.pass_3d = (rmsd is not None and rmsd < rmsd_thresh)
    rmsd_str = f"{rmsd:.2f} Å" if rmsd is not None else "N/A"
    print(f"✔ Boltz-2 | RMSD={rmsd_str} {'✔' if probe.pass_3d else '✘'}")
    return probe

# ─── 9. Outputs por gene ──────────────────────────────────────────────────────
def write_gene_outputs(probes: list[Probe], gene_key: str):
    out_dir = ALIGN_DIR / gene_key
    tsv     = out_dir / f"{gene_key}_probes_scored.tsv"
    fasta   = out_dir / f"{gene_key}_viroscope_probes.fasta"
    fields  = list(asdict(probes[0]).keys())
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

# ─── 10. CSV consolidado (Romeu + Beatriz) ────────────────────────────────────
def _parse_beatriz_xlsx(path: Path) -> list[dict]:
    """
    Lê o XLSX da Beatriz e normaliza para o schema do pipeline.
    Sheet 1 "Nossas sequencias": cabeçalho na linha 2 (0-indexed), dados a partir da linha 3.
      Colunas-chave: Column1, species, Sequence, Gene/STRAIN, Tm, GC content,
                     Lowest folding energy kcal/mole at 25C, notes
    Sheet 2 "New probes (papers)": cabeçalho na linha 0.
      Colunas-chave: species, Published probes (5'-3'), gene, notes, Source
    Sequências de literatura: Tm/GC/hairpin calculados via primer3 in-situ.
    """
    rows: list[dict] = []
    if not path.exists():
        print(f"  ⚠ Ficheiro Beatriz não encontrado: {path}")
        return rows

    _NORM_ORG = {
        "Strepptocuccus": "Streptococcus",
        "Strepptococcus": "Streptococcus",
    }
    _GROUP_A  = ["Staphylococcus aureus", "Neisseria meningitidis"]

    def _clean_seq(s) -> str:
        if pd.isna(s):
            return ""
        s = str(s).upper()
        s = re.sub(r"^[A-Z0-9]+-C\d+-", "", s)   # strip NH2-C6-, etc.
        return re.sub(r"[^ATGCUN]", "", s)

    def _parse_float_unit(s) -> Optional[float]:
        if pd.isna(s):
            return None
        m = re.search(r"[-\d.]+", str(s))
        return float(m.group()) if m else None

    def _parse_gc_pct(s) -> Optional[float]:
        v = _parse_float_unit(s)
        return round(v / 100, 3) if v is not None else None

    def _norm_org(s) -> str:
        if pd.isna(s):
            return ""
        s = str(s).strip()
        for wrong, right in _NORM_ORG.items():
            s = s.replace(wrong, right)
        return s

    def _extract_gene(s) -> str:
        if pd.isna(s):
            return ""
        s = str(s).strip()
        return s.split(":")[0].strip() if ":" in s else s.split()[0]

    def _group(org: str) -> str:
        return "A" if any(g in org for g in _GROUP_A) else "B"

    def _thresh(gene: str, param: str):
        gk = gene.lower() if gene.lower() in TARGETS else None
        return TARGETS[gk].get(param, DEFAULTS[param]) if gk else DEFAULTS[param]

    def _pass_basic(tm, gc, hp, gene: str) -> Optional[bool]:
        if tm is None or gc is None:
            return None
        return bool(
            _thresh(gene, "tm_min") <= tm <= _thresh(gene, "tm_max") and
            _thresh(gene, "gc_min") <= gc <= _thresh(gene, "gc_max") and
            (hp is None or hp >= _thresh(gene, "hp_min"))
        )

    xl = pd.ExcelFile(path)

    # ── Sheet 1: "Nossas sequencias" ─────────────────────────────────────────
    # Detecta o cabeçalho dinamicamente (linha que contém "Sequence")
    _raw1 = pd.read_excel(xl, sheet_name="Nossas sequencias", header=None)
    _hrow = next(
        (i for i, r in _raw1.iterrows()
         if any(str(v).strip().lower() == "sequence" for v in r.values if pd.notna(v))),
        2
    )
    df1 = pd.read_excel(xl, sheet_name="Nossas sequencias",
                        header=_hrow, dtype=str)
    df1.columns = [str(c).strip() for c in df1.columns]

    for _, row in df1.iterrows():
        seq = _clean_seq(row.get("Sequence"))
        if len(seq) < 10:
            continue
        organism = _norm_org(row.get("species"))
        gene     = _extract_gene(row.get("Gene/STRAIN"))
        tm       = _parse_float_unit(row.get("Tm"))
        gc       = _parse_gc_pct(row.get("GC content"))
        hp       = _parse_float_unit(row.get("Lowest folding energy kcal/mole at 25C"))
        pid_raw  = row.get("Column1")
        pid      = str(pid_raw).strip() if pd.notna(pid_raw) else f"bea{len(rows)+1:03d}"
        rows.append({
            "probe_id":      pid,
            "gene":          gene,
            "organism":      organism,
            "group":         _group(organism),
            "sequence":      seq,
            "pos_start":     "",
            "pos_end":       "",
            "tm":            tm if tm is not None else "",
            "gc":            gc if gc is not None else "",
            "hairpin_dg":    hp,
            "homodimer_dg":  "",
            "nupack_dg":     "",
            "nupack_defect": "",
            "structure_cif": "",
            "rmsd":          "",
            "pass_basic":    _pass_basic(tm, gc, hp, gene),
            "pass_nupack":   "",
            "pass_3d":       "",
            "fonte":         "beatriz",
            "notes":         str(row.get("notes", "") or "").strip(),
        })

    # ── Sheet 2: "New probes (papers)" ───────────────────────────────────────
    _raw2 = pd.read_excel(xl, sheet_name="New probes (papers)", header=None)
    _hrow2 = next(
        (i for i, r in _raw2.iterrows()
         if any("probe" in str(v).lower() or "species" in str(v).lower()
                for v in r.values if pd.notna(v))),
        0
    )
    df2 = pd.read_excel(xl, sheet_name="New probes (papers)",
                        header=_hrow2, dtype=str)
    df2.columns = [str(c).strip() for c in df2.columns]
    for _, row in df2.iterrows():
        raw  = str(row.get("Published probes (5'-3')", "") or "")
        seq  = re.sub(r"[^ATGCUN]", "", raw.upper())
        if len(seq) < 10:
            continue
        organism = _norm_org(row.get("species"))
        gene     = str(row.get("gene", "") or "").strip()
        source   = str(row.get("Source", "") or "").strip()
        notes    = str(row.get("notes", "") or "").strip()
        if source:
            notes = f"{notes} [Source: {source}]".strip("[ ]")
        sc   = score_probe(seq)
        pid  = (f"lit{len(rows)+1:03d}_"
                f"{_ORG_ABBR.get(organism, organism[:4].capitalize())}_"
                f"{gene}")
        rows.append({
            "probe_id":      pid,
            "gene":          gene,
            "organism":      organism,
            "group":         _group(organism),
            "sequence":      seq,
            "pos_start":     "",
            "pos_end":       "",
            "tm":            sc["tm"],
            "gc":            sc["gc"],
            "hairpin_dg":    sc["hairpin_dg"],
            "homodimer_dg":  sc["homodimer_dg"],
            "nupack_dg":     "",
            "nupack_defect": "",
            "structure_cif": "",
            "rmsd":          "",
            "pass_basic":    _pass_basic(sc["tm"], sc["gc"], sc["hairpin_dg"], gene),
            "pass_nupack":   "",
            "pass_3d":       "",
            "fonte":         "literatura",
            "notes":         notes,
        })

    print(f"  ✔ Beatriz: {len(rows)} probes lidas do XLSX")
    return rows

def write_consolidated_csv(all_probes: list[Probe]):
    """
    Gera FINAL_PROBES_ALL.csv com probes do Romeu + Beatriz,
    ordenado por gene e por pass_basic desc.
    """
    out_path = BASE_DIR / "output" / "FINAL_PROBES_ALL.csv"
    romeu_rows   = [asdict(p) for p in all_probes]
    beatriz_rows = _parse_beatriz_xlsx(BEATRIZ_XLSX)
    all_rows     = romeu_rows + beatriz_rows
    all_rows.sort(key=lambda r: (r["gene"], not r["pass_basic"]))
    fields = list(asdict(Probe("","","")).keys())
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"\n  ✔ CSV consolidado: {out_path}  ({len(all_rows)} probes total)")

# ─── Pipeline principal ───────────────────────────────────────────────────────
def run_pipeline(run_nupack=True, run_3d=True):
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

        # 1. Download
        records = fetch_sequences(gene_key)

        # 2. Alinhamento
        aln = align_mafft(records, gene_key)

        # 3. Janelas
        windows = candidate_windows(aln, gene_key)

        # 4. Scoring básico
        probes: list[Probe] = []
        print(f"  [4/6] Scoring básico (Tm, GC, hairpin, homodimer)...")
        for w in windows:
            sc = score_probe(w["seq"])
            p  = Probe(
                gene       = gene_key,
                organism   = t["organism"],
                group      = t["group"],
                sequence   = w["seq"],
                pos_start  = w["start"],
                pos_end    = w["end"],
                tm         = sc["tm"],
                gc         = sc["gc"],
                hairpin_dg = sc["hairpin_dg"],
                homodimer_dg = sc["homodimer_dg"],
            )
            p.pass_basic = passes_basic(p, gene_key)
            p.probe_id   = build_probe_id(p, next(probe_counter))
            probes.append(p)

        n_pass = sum(p.pass_basic for p in probes)
        print(f"    ✔ {n_pass}/{len(probes)} probes passam critérios básicos")

        # 5. seqfold — estrutura secundária (só nas probes básicas aprovadas)
        if run_nupack:
            print(f"  [5/6] seqfold (ΔG MFE + fracção emparelhada) ...")
            for p in probes:
                if p.pass_basic:
                    p = run_seqfold_probe(p, gene_key)
            n_nup = sum(p.pass_nupack for p in probes if p.pass_basic)
            print(f"    ✔ {n_nup}/{n_pass} probes passam seqfold")

        # 6. Estrutura 3D + RMSD (só PASS básico + seqfold)
        if run_3d:
            _nf_ok     = NUCLEOFOLD3D_SCRIPT.exists()
            _boltz_ok  = _boltz_cmd() is not None
            if not _nf_ok and not _boltz_ok:
                print(f"  [6/6] Estrutura 3D: ferramentas não disponíveis — a saltar.")
                print(f"    NucleoFold3D : não encontrado em {NUCLEOFOLD3D_SCRIPT}")
                print(f"    Boltz-2      : não instalado neste Python ({sys.executable})")
                print(f"    → Para 3D: pip install boltz  ou  correr com python conda")
            else:
                print(f"  [6/6] Estrutura 3D (NucleoFold → Boltz-2)...")
                for p in probes:
                    gate = p.pass_basic and (p.pass_nupack if run_nupack else True)
                    if gate:
                        p = run_3d_pipeline(p)

        write_gene_outputs(probes, gene_key)
        all_probes.extend(probes)

    # CSV final consolidado
    write_consolidated_csv(all_probes)

    # Resumo
    print("\n" + "═"*60)
    print("  RESUMO FINAL")
    print("═"*60)
    for gene_key in TARGETS:
        ps = [p for p in all_probes if p.gene == gene_key]
        nb = sum(p.pass_basic for p in ps)
        nn = sum(p.pass_nupack for p in ps)
        n3 = sum(p.pass_3d for p in ps)
        print(f"  {gene_key:<8} basic={nb}  nupack={nn}  3d={n3}")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="GFET Probe Pipeline v4")
    ap.add_argument("--no-nupack", action="store_true", help="Saltar seqfold (estrutura secundária)")
    ap.add_argument("--no-3d",    action="store_true", help="Saltar modelação 3D")
    args = ap.parse_args()
    run_pipeline(run_nupack=not args.no_nupack, run_3d=not args.no_3d)
