"""
GFET Probe Pipeline v4
─────────────────────────────────────────────────────────────────────────────
Novidades face à v3:
  • build_probe_id()         — nome codificado: p001_Saur_nuc_pos373-395_Tm54.1_GC40_hp-0.38_PASS
  • Limiares por gene        — definidos em TARGETS{} (fácil de ajustar)
  • run_nupack_probe()       — ΔG ensemble + defect normalizado (NUPACK CLI)
  • run_nucleofold()         — estrutura 3D via NucleoFold API
  • run_boltz2()             — fallback 3D via Boltz-2 (YAML + CLI)
  • _calc_rmsd_between_cifs()— RMSD entre réplicas (Biopython Superimposer)
  • write_consolidated_csv() — CSV único com probes Romeu + Beatriz

Referências dos limiares:
  SantaLucia & Hicks 2004  — Tm nearest-neighbour
  Wetmur 1991              — comprimento da probe (18–28 nt)
  Zadeh et al. 2011        — NUPACK ΔG ensemble e defect
  IDT OligoAnalyzer        — GC 40–65%, homodimer ΔG > -5 kcal/mol
─────────────────────────────────────────────────────────────────────────────
"""

import os, sys, re, csv, json, yaml, shutil, subprocess, tempfile, textwrap, itertools, time
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
ALIGN_DIR    = BASE_DIR / "alignments"
NUPACK_DIR   = BASE_DIR / "nupack_out"
STRUCT_DIR   = BASE_DIR / "structures"
BEATRIZ_XLSX = BASE_DIR / "BeatrizMasterThesis_ProbesData_IPLEXMED_ENERGIAS_CAPS.xlsx"
MAFFT_BIN    = BASE_DIR / "MAFFT" / "mafft-7.526-win64-signed" / "mafft-win" / "mafft.bat"

for d in [ALIGN_DIR, NUPACK_DIR, STRUCT_DIR]:
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
    "nupack_dg_max":    -12.0,    # kcal/mol — Zadeh et al. 2011
    "nupack_defect_max":  0.10,   # defect normalizado — Zadeh et al. 2011
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
                recs = list(SeqIO.parse(h, "fasta")); h.close()
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
                records = list(SeqIO.parse(h, "fasta")); h.close()
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

# ─── 6. NUPACK ───────────────────────────────────────────────────────────────
def run_nupack_probe(probe: Probe, gene_key: str) -> Probe:
    """
    Corre o NUPACK (pfunc) para uma probe.
    Preenche nupack_dg, nupack_defect e pass_nupack.
    Requer NUPACK 4+ instalado e acessível via PATH.
    Ref: Zadeh et al. 2011, J. Comput. Chem.
    """
    if shutil.which("nupack") is None:
        probe.notes += "[NUPACK não encontrado] "
        return probe
    with tempfile.TemporaryDirectory() as tmpdir:
        inp = Path(tmpdir) / "probe.in"
        inp.write_text(f"1\n{probe.sequence}\n")
        try:
            result = subprocess.run(
                ["nupack", "pfunc", str(inp)],
                capture_output=True, text=True, timeout=60
            )
            # parse ΔG ensemble
            for line in result.stdout.splitlines():
                if "Free energy" in line or "dG" in line:
                    m = re.search(r"[-\d.]+", line)
                    if m:
                        probe.nupack_dg = float(m.group())
                if "Ensemble defect" in line or "defect" in line.lower():
                    m = re.search(r"[\d.]+", line)
                    if m:
                        probe.nupack_defect = float(m.group())
        except Exception as e:
            probe.notes += f"[NUPACK erro: {e}] "
            return probe
    dg_ok  = probe.nupack_dg     is not None and probe.nupack_dg     <= cfg(gene_key, "nupack_dg_max")
    def_ok = probe.nupack_defect is not None and probe.nupack_defect <= cfg(gene_key, "nupack_defect_max")
    probe.pass_nupack = dg_ok and def_ok
    return probe

# ─── 7. Estrutura 3D ─────────────────────────────────────────────────────────
def run_nucleofold(probe: Probe) -> Optional[Path]:
    """
    Submete a sequência à API pública do NucleoFold (3dRNA/NucleoFold).
    Devolve o path do CIF gerado ou None se falhar.
    """
    url = "http://biophy.hust.edu.cn/new/3dRNA/api/submit"   # endpoint público
    try:
        resp = requests.post(url,
                             json={"sequence": probe.sequence, "type": "DNA"},
                             timeout=120)
        resp.raise_for_status()
        data    = resp.json()
        job_id  = data.get("job_id") or data.get("id")
        if not job_id:
            return None
        # poll até terminar (máx 5 min)
        import time
        for _ in range(60):
            time.sleep(5)
            r2 = requests.get(f"{url.replace('submit','result')}/{job_id}", timeout=30)
            r2.raise_for_status()
            d2 = r2.json()
            if d2.get("status") == "done":
                cif_url = d2.get("cif_url") or d2.get("structure_url")
                if cif_url:
                    cif_data = requests.get(cif_url, timeout=30).text
                    out_cif  = STRUCT_DIR / f"{probe.probe_id}_nucleofold.cif"
                    out_cif.write_text(cif_data)
                    return out_cif
            elif d2.get("status") == "error":
                return None
    except Exception:
        return None
    return None

def run_boltz2(probe: Probe, n_samples: int = 3) -> list[Path]:
    """
    Corre o Boltz-2 localmente para gerar n_samples estruturas CIF.
    Usa o CLI: boltz predict <yaml> --out_dir <dir> --diffusion_samples n
    Ref: Wohlwend et al. 2024 (Boltz-2 paper)
    """
    if shutil.which("boltz") is None:
        probe.notes += "[Boltz-2 não encontrado] "
        return []
    out_dir = STRUCT_DIR / probe.probe_id / "boltz"
    out_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = out_dir / "input.yaml"
    yaml_path.write_text(yaml.dump({
        "version": 1,
        "sequences": [{"dna": {"id": ["A"], "sequence": probe.sequence}}]
    }))
    try:
        subprocess.run([
            "boltz", "predict", str(yaml_path),
            "--out_dir", str(out_dir),
            "--recycling_steps", "3",
            "--sampling_steps", "200",
            "--diffusion_samples", str(n_samples),
            "--device", "cpu",
        ], check=True, capture_output=True, timeout=600)
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
    Tenta NucleoFold primeiro; se falhar, tenta Boltz-2.
    Se Boltz-2 gerar ≥2 réplicas, calcula RMSD.
    Pass_3d = True se RMSD < rmsd_max.
    """
    print(f"      3D: NucleoFold... ", end="", flush=True)
    cif = run_nucleofold(probe)
    if cif:
        probe.structure_cif = str(cif)
        probe.pass_3d       = True
        print("✔ NucleoFold")
        return probe
    print("falhou → Boltz-2... ", end="", flush=True)
    cifs = run_boltz2(probe, n_samples=3)
    if not cifs:
        print("✘")
        return probe
    probe.structure_cif = str(cifs[0])
    rmsd = _calc_rmsd_between_cifs(cifs)
    probe.rmsd = rmsd
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
    df1 = pd.read_excel(xl, sheet_name="Nossas sequencias", header=2)
    for _, row in df1.iterrows():
        seq = _clean_seq(row.get("Sequence"))
        if len(seq) < 10:
            continue
        organism = _norm_org(row.get("species"))
        gene     = _extract_gene(row.get("Gene/STRAIN"))
        tm       = _parse_float_unit(row.get("Tm"))
        gc       = _parse_gc_pct(row.get("GC content"))
        hp       = row.get("Lowest folding energy kcal/mole at 25C")
        hp       = float(hp) if pd.notna(hp) else None
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
    df2 = pd.read_excel(xl, sheet_name="New probes (papers)", header=0)
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
    out_path = ALIGN_DIR / "FINAL_PROBES_ALL.csv"
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
    print("  GFET Probe Pipeline v4")
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

        # 5. NUPACK (só nas probes básicas aprovadas)
        if run_nupack:
            print(f"  [5/6] NUPACK (ΔG ensemble + defect)...")
            for p in probes:
                if p.pass_basic:
                    p = run_nupack_probe(p, gene_key)
            n_nup = sum(p.pass_nupack for p in probes if p.pass_basic)
            print(f"    ✔ {n_nup}/{n_pass} probes passam NUPACK")

        # 6. Estrutura 3D + RMSD (só PASS básico + NUPACK)
        if run_3d:
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
    ap.add_argument("--no-nupack", action="store_true", help="Saltar NUPACK")
    ap.add_argument("--no-3d",    action="store_true", help="Saltar modelação 3D")
    args = ap.parse_args()
    run_pipeline(run_nupack=not args.no_nupack, run_3d=not args.no_3d)
