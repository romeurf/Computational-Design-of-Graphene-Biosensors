#!/usr/bin/env python3
"""
analysis.py — Análise exploratória / estatística do pipeline GFET.

Um único ponto de entrada para toda a análise pedida pelo orientador:
  1. seqfold  — distribuição do ΔG MFE + curva cumulativa de aprovação ao longo de
                TODO o intervalo (para definir um limiar de forma fundamentada).
  2. tamanhos — nº de sequências NCBI (hits) e tamanho da região/gene por gene,
                comprimento das probes, com teste de normalidade (Shapiro-Wilk).
  3. diversidade — grau de diversidade entre as sequências recuperadas SEM
                   alinhamento (distância k-mer par-a-par).
  4. rarefação — diversidade/riqueza de k-mers vs nº de sequências → justifica o N
                 a usar (escalar com rigor, não "até ao infinito").
  5. descritores — alternativa ao PyBioMed: descritores de sequência por probe.

Saídas:
  output/figures/*.png           (figuras)
  output/analysis/analysis_summary.md   (resumo legível)
  output/analysis/seqfold_threshold_sweep.csv
  output/analysis/per_gene_sizes.csv
  output/analysis/probe_descriptors.csv

Uso:  python scripts/analysis.py            (corre tudo)
"""

import csv
import itertools
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics.pairwise import cosine_distances

# ── Paths / config ────────────────────────────────────────────────────────────
BASE   = Path(__file__).resolve().parent.parent
ALIGN  = BASE / "output" / "alignments"
FINAL  = BASE / "output" / "FINAL_PROBES_ALL.csv"
FIGDIR = BASE / "output" / "figures"
ANADIR = BASE / "output" / "analysis"
FIGDIR.mkdir(parents=True, exist_ok=True)
ANADIR.mkdir(parents=True, exist_ok=True)

GENES = ["nuc", "rmpM", "lytA", "oprL", "algD", "frdB"]
SEQFOLD_DEFAULT = -6.0           # limiar atual no pipeline (Zadeh et al. 2011)
KMER_K = 4

_md_lines: list[str] = []
def md(line: str = ""):
    _md_lines.append(line)
    print(line)

# ── Leitura ───────────────────────────────────────────────────────────────────
def read_fasta(path: Path) -> list[str]:
    if not path.exists():
        return []
    seqs, buf, have = [], [], False
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(">"):
            if have:
                seqs.append("".join(buf))
            buf, have = [], True
        elif have:
            buf.append(line.strip())
    if have:
        seqs.append("".join(buf))
    return [s.upper() for s in seqs if s.strip()]

def gene_input_seqs(gene: str) -> list[str]:
    """Sequências NCBI cruas; fallback p/ aligned.fasta (sem gaps) se input ausente."""
    seqs = read_fasta(ALIGN / gene / "input.fasta")
    if not seqs:
        seqs = [s.replace("-", "") for s in read_fasta(ALIGN / gene / "aligned.fasta")]
    return [s for s in seqs if s]

def load_probes() -> pd.DataFrame:
    df = pd.read_csv(FINAL)
    for c in ("tm", "gc", "conservation", "seqfold_dg", "seqfold_paired_frac",
              "hairpin_dg", "homodimer_dg", "nofold_score"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("pass_basic", "pass_seqfold"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.lower().isin(("true", "1", "1.0"))
    return df

# ── 1. seqfold: distribuição + curva cumulativa de aprovação ───────────────────
def analyze_seqfold(df: pd.DataFrame):
    md("\n## 1. seqfold — distribuição de ΔG MFE e limiar")
    # População relevante: probes próprias que passaram o filtro básico e têm ΔG
    pop = df[(df["fonte"] == "romeu") & (df["pass_basic"]) & (df["seqfold_dg"].notna())]
    dg = pop["seqfold_dg"].to_numpy()
    if dg.size == 0:
        md("_Sem dados de seqfold._"); return
    md(f"- Probes consideradas (próprias, pass_basic, ΔG definido): **{dg.size}**")
    md(f"- ΔG MFE: min **{dg.min():.2f}**, mediana **{np.median(dg):.2f}**, "
       f"média **{dg.mean():.2f}**, máx **{dg.max():.2f}** kcal/mol")
    pcts = {p: float(np.percentile(dg, p)) for p in (1, 5, 10, 25, 50)}
    md("- Percentis de ΔG (cauda mais estável = mais negativa): "
       + ", ".join(f"P{p}={v:.2f}" for p, v in pcts.items()))

    # Histograma global + por gene
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    axes[0].hist(dg, bins=40, color="#2E75B6", edgecolor="white")
    axes[0].axvline(SEQFOLD_DEFAULT, color="red", ls="--",
                    label=f"limiar atual {SEQFOLD_DEFAULT}")
    axes[0].set(title="Distribuição ΔG MFE (todas as probes pass_basic)",
                xlabel="ΔG (kcal/mol)", ylabel="nº de probes")
    axes[0].legend()
    for g in GENES:
        gd = pop[pop["gene"] == g]["seqfold_dg"].to_numpy()
        if gd.size > 5:
            axes[1].hist(gd, bins=25, histtype="step", lw=1.6, label=f"{g} (n={gd.size})")
    axes[1].axvline(SEQFOLD_DEFAULT, color="red", ls="--")
    axes[1].set(title="ΔG MFE por gene", xlabel="ΔG (kcal/mol)", ylabel="nº de probes")
    axes[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIGDIR / "seqfold_hist.png", dpi=130); plt.close(fig)

    # Curva cumulativa de aprovação ao longo de TODO o intervalo de ΔG
    grid = np.linspace(dg.min() - 0.5, dg.max() + 0.5, 300)
    frac_pass = np.array([(dg >= t).mean() for t in grid])  # ΔG >= t → passa
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.plot(grid, frac_pass * 100, color="#1F3864", lw=2)
    ax.axvline(SEQFOLD_DEFAULT, color="red", ls="--",
               label=f"limiar atual {SEQFOLD_DEFAULT} → {(dg>=SEQFOLD_DEFAULT).mean()*100:.1f}% passam")
    for p in (5, 10, 25):
        v = np.percentile(dg, p)
        ax.axvline(v, color="grey", ls=":", lw=0.9)
        ax.text(v, 5, f"P{p}", rotation=90, fontsize=7, color="grey", va="bottom")
    ax.set(title="Curva cumulativa de aprovação seqfold (% que passa para cada limiar)",
           xlabel="limiar de ΔG (kcal/mol)  — probe passa se ΔG ≥ limiar",
           ylabel="% de probes que passam")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(FIGDIR / "seqfold_cumulative_pass.png", dpi=130); plt.close(fig)

    # Tabela de sweep de limiares (global + por gene)
    thresholds = [0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0]
    rows = []
    for t in thresholds:
        row = {"threshold": t, "global_pass": int((dg >= t).sum()),
               "global_pct": round((dg >= t).mean() * 100, 1)}
        for g in GENES:
            gd = pop[pop["gene"] == g]["seqfold_dg"].to_numpy()
            row[f"{g}_pass"] = int((gd >= t).sum()) if gd.size else 0
        rows.append(row)
    pd.DataFrame(rows).to_csv(ANADIR / "seqfold_threshold_sweep.csv", index=False)
    md("\n- Aprovação a vários limiares (ΔG ≥ limiar):")
    md("\n| limiar | passam (global) | % | " + " | ".join(GENES) + " |")
    md("|---|---|---|" + "---|" * len(GENES))
    for r in rows:
        md(f"| {r['threshold']:.0f} | {r['global_pass']} | {r['global_pct']} | "
           + " | ".join(str(r[f'{g}_pass']) for g in GENES) + " |")

    # Recomendação fundamentada
    rec = float(np.percentile(dg, 5))
    md(f"\n**Recomendação (para discussão):** o limiar atual {SEQFOLD_DEFAULT} kcal/mol é "
       f"muito permissivo — deixa passar {(dg>=SEQFOLD_DEFAULT).mean()*100:.1f}% das probes. "
       f"Para remover apenas as ~5% com estrutura secundária mais estável (cauda mais negativa), "
       f"um limiar ≈ **{rec:.1f} kcal/mol** (P5 da distribuição) é defensável estatisticamente. "
       f"NÃO foi alterado no pipeline — é uma recomendação baseada na distribuição observada.")

# ── 2. tamanhos & hits por gene + normalidade ──────────────────────────────────
def analyze_sizes(df: pd.DataFrame):
    md("\n## 2. Tamanhos & hits por gene + normalidade")
    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    for ax, g in zip(axes.ravel(), GENES):
        seqs = gene_input_seqs(g)
        lens = np.array([len(s) for s in seqs])
        n = len(lens)
        if n == 0:
            ax.set_title(f"{g} — sem sequências"); ax.axis("off")
            rows.append({"gene": g, "n_seqs": 0}); continue
        sw_p = float(stats.shapiro(lens).pvalue) if 3 <= n <= 5000 and lens.std() > 0 else float("nan")
        normal = (sw_p >= 0.05) if not math.isnan(sw_p) else None
        rows.append({"gene": g, "n_seqs": n,
                     "len_min": int(lens.min()), "len_max": int(lens.max()),
                     "len_mean": round(float(lens.mean()), 1),
                     "len_sd": round(float(lens.std(ddof=1)) if n > 1 else 0.0, 1),
                     "shapiro_p": round(sw_p, 4) if not math.isnan(sw_p) else "n/a",
                     "normal_p<0.05": ("sim" if normal else "não") if normal is not None else "n/a"})
        ax.hist(lens, bins=min(30, max(5, n // 3)), color="#70AD47", edgecolor="white")
        ax.set_title(f"{g}  (n={n}, μ={lens.mean():.0f}±{lens.std():.0f})", fontsize=10)
        ax.set_xlabel("comprimento (bp)")
    fig.suptitle("Distribuição do tamanho das sequências NCBI recuperadas por gene")
    fig.tight_layout(); fig.savefig(FIGDIR / "sizes_per_gene.png", dpi=130); plt.close(fig)

    pd.DataFrame(rows).to_csv(ANADIR / "per_gene_sizes.csv", index=False)
    md("\n| gene | nº seqs (hits) | comp. min–máx | média±DP | Shapiro p | normal? |")
    md("|---|---|---|---|---|---|")
    for r in rows:
        if r.get("n_seqs", 0) == 0:
            md(f"| {r['gene']} | 0 | — | — | — | — |"); continue
        md(f"| {r['gene']} | {r['n_seqs']} | {r['len_min']}–{r['len_max']} | "
           f"{r['len_mean']}±{r['len_sd']} | {r['shapiro_p']} | {r['normal_p<0.05']} |")
    md("\n_Shapiro-Wilk: p < 0.05 ⇒ rejeita normalidade. n/a quando n<3 ou variância nula._")

    # QQ-plots dos tamanhos por gene (uma figura)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    for ax, g in zip(axes.ravel(), GENES):
        lens = np.array([len(s) for s in gene_input_seqs(g)])
        if lens.size >= 3 and lens.std() > 0:
            stats.probplot(lens, dist="norm", plot=ax)
            ax.set_title(f"QQ — {g}", fontsize=10)
        else:
            ax.set_title(f"{g} — n/a"); ax.axis("off")
    fig.suptitle("QQ-plots (normalidade do tamanho das sequências por gene)")
    fig.tight_layout(); fig.savefig(FIGDIR / "sizes_qqplots.png", dpi=130); plt.close(fig)

    # Comprimento das probes (próprias)
    plen = df[df["fonte"] == "romeu"]["sequence"].dropna().map(len).to_numpy()
    if plen.size:
        md(f"\n- Comprimento das probes próprias: n={plen.size}, "
           f"intervalo {plen.min()}–{plen.max()} nt, média {plen.mean():.1f} "
           f"(discreto/limitado 18–28 nt → não-normal por construção).")
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(plen, bins=range(int(plen.min()), int(plen.max()) + 2),
                color="#2E75B6", edgecolor="white", align="left")
        ax.set(title="Comprimento das probes (nt)", xlabel="nt", ylabel="nº de probes")
        fig.tight_layout(); fig.savefig(FIGDIR / "probe_length.png", dpi=130); plt.close(fig)

# ── 3. diversidade sem alinhamento (k-mer) ──────────────────────────────────────
def kmer_matrix(seqs: list[str], k: int = KMER_K):
    kmers = ["".join(p) for p in itertools.product("ACGT", repeat=k)]
    idx = {km: i for i, km in enumerate(kmers)}
    M = np.zeros((len(seqs), len(kmers)), dtype=float)
    for r, s in enumerate(seqs):
        s = "".join(c for c in s if c in "ACGT")
        tot = 0
        for i in range(len(s) - k + 1):
            j = idx.get(s[i:i + k])
            if j is not None:
                M[r, j] += 1; tot += 1
        if tot:
            M[r] /= tot
    return M

def analyze_diversity(df: pd.DataFrame):
    md("\n## 3. Diversidade entre sequências recuperadas (sem alinhamento, k-mer)")
    md(f"- Distância = 1 − similaridade do cosseno entre vetores de frequência de {KMER_K}-mers.")
    rows, heat = [], {}
    for g in GENES:
        seqs = gene_input_seqs(g)
        if len(seqs) < 2:
            rows.append((g, len(seqs), "n/a")); continue
        M = kmer_matrix(seqs)
        D = cosine_distances(M)
        mean_div = float(D[np.triu_indices(len(seqs), 1)].mean())
        rows.append((g, len(seqs), round(mean_div, 4)))
        heat[g] = D
    md("\n| gene | nº seqs | diversidade média (0=idênticas, →1 diversas) |")
    md("|---|---|---|")
    for g, n, d in rows:
        md(f"| {g} | {n} | {d} |")

    # Heatmaps (subamostra até 40 seqs por gene para legibilidade)
    valid = [g for g in GENES if g in heat]
    if valid:
        fig, axes = plt.subplots(1, len(valid), figsize=(3.2 * len(valid), 3.4))
        if len(valid) == 1:
            axes = [axes]
        for ax, g in zip(axes, valid):
            D = heat[g]
            if D.shape[0] > 40:
                sel = np.linspace(0, D.shape[0] - 1, 40).astype(int)
                D = D[np.ix_(sel, sel)]
            im = ax.imshow(D, cmap="viridis", vmin=0)
            ax.set_title(f"{g}", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=axes, fraction=0.025, label="distância k-mer")
        fig.suptitle("Diversidade par-a-par (distância k-mer, sem alinhamento)")
        fig.savefig(FIGDIR / "diversity_heatmaps.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ── 4. rarefação → justificar N ─────────────────────────────────────────────────
def analyze_rarefaction(df: pd.DataFrame):
    md("\n## 4. Rarefação — quantas sequências usar (escalar com rigor)")
    md("- Para cada gene, subamostram-se N sequências e mede-se a diversidade k-mer média "
       "e a riqueza (k-mers distintos). Quando saturam, mais sequências acrescentam pouco.")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    rng = np.random.default_rng(0)
    sat_rows = []
    for g in GENES:
        seqs = gene_input_seqs(g)
        n = len(seqs)
        if n < 5:
            continue
        M = kmer_matrix(seqs)
        sched = [m for m in (5, 10, 20, 30, 50, 75, 100, 150, 200) if m <= n]
        if n not in sched:
            sched.append(n)
        divs, richs = [], []
        for m in sched:
            dd, rr = [], []
            for _ in range(5):
                sel = rng.choice(n, m, replace=False)
                sub = M[sel]
                if m >= 2:
                    D = cosine_distances(sub)
                    dd.append(D[np.triu_indices(m, 1)].mean())
                rr.append(int((sub.sum(axis=0) > 0).sum()))
            divs.append(np.mean(dd) if dd else 0.0)
            richs.append(np.mean(rr))
        axes[0].plot(sched, divs, marker="o", ms=3, label=g)
        axes[1].plot(sched, richs, marker="o", ms=3, label=g)
        # N de saturação: primeiro m onde riqueza ≥ 95% da riqueza máxima
        rich_max = richs[-1]
        sat_n = next((m for m, r in zip(sched, richs) if rich_max and r >= 0.95 * rich_max), sched[-1])
        sat_rows.append((g, n, sat_n))
    axes[0].set(title="Diversidade k-mer média vs N", xlabel="nº de sequências", ylabel="diversidade média")
    axes[1].set(title="Riqueza (k-mers distintos) vs N", xlabel="nº de sequências", ylabel="k-mers distintos")
    for ax in axes:
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(FIGDIR / "rarefaction.png", dpi=130); plt.close(fig)

    md("\n| gene | nº disponível | N de saturação (riqueza ≥95% do máx) |")
    md("|---|---|---|")
    for g, n, s in sat_rows:
        md(f"| {g} | {n} | {s} |")
    if sat_rows:
        rec_n = max(s for _, _, s in sat_rows)
        md(f"\n**Recomendação de N:** a diversidade/riqueza satura por volta de "
           f"**N ≈ {rec_n}** sequências para o gene mais exigente. Usar N ≈ {rec_n}–{max(rec_n,100)} "
           f"por gene é suficiente (mais do que isso acrescenta pouca informação nova). "
           f"Genes com poucas sequências no NCBI (ex.: algD) são o fator limitante real, não o cap.")

# ── 5. descritores de sequência (alternativa ao PyBioMed) ───────────────────────
def seq_descriptors(seq: str) -> dict:
    s = "".join(c for c in seq.upper() if c in "ACGT")
    L = len(s) or 1
    a, c, g, t = (s.count(x) / L for x in "ACGT")
    pur = (s.count("A") + s.count("G")) / L
    ent = -sum(p * math.log2(p) for p in (a, c, g, t) if p > 0)
    d = {"length": len(s), "pA": round(a, 3), "pC": round(c, 3), "pG": round(g, 3),
         "pT": round(t, 3), "GC": round(g + c, 3), "AT": round(a + t, 3),
         "purine": round(pur, 3), "pyrimidine": round(1 - pur, 3), "entropy": round(ent, 3)}
    dinucs = ["".join(p) for p in itertools.product("ACGT", repeat=2)]
    cnt = {dn: 0 for dn in dinucs}
    for i in range(len(s) - 1):
        cnt[s[i:i + 2]] += 1
    tot = max(len(s) - 1, 1)
    for dn in dinucs:
        d[f"di_{dn}"] = round(cnt[dn] / tot, 4)
    return d

def compute_descriptors(df: pd.DataFrame):
    md("\n## 5. Descritores de sequência (alternativa ao PyBioMed)")
    sub = df[df["sequence"].notna()].copy()
    recs = []
    for _, row in sub.iterrows():
        base = {"probe_id": row.get("probe_id", ""), "gene": row.get("gene", ""),
                "fonte": row.get("fonte", "")}
        base.update(seq_descriptors(str(row["sequence"])))
        recs.append(base)
    out = pd.DataFrame(recs)
    out.to_csv(ANADIR / "probe_descriptors.csv", index=False)
    md(f"- {len(out)} probes descritas → `output/analysis/probe_descriptors.csv` "
       f"({out.shape[1]} colunas: composição, GC/AT, purina/pirimidina, entropia, "
       f"16 dinucleótidos).")
    md("- _Descritores de estruturas 3D: adiados (sem estruturas em disco) — ficam como próximo passo._")

# ── main ────────────────────────────────────────────────────────────────────────
def main():
    if not FINAL.exists():
        raise SystemExit(f"Não encontrado: {FINAL} — corre primeiro o pipeline.")
    df = load_probes()
    md("# Análise exploratória — pipeline GFET")
    md(f"Fonte: `{FINAL.name}` ({len(df)} probes; "
       f"{int((df['fonte']=='romeu').sum())} próprias, "
       f"{int((df['fonte']!='romeu').sum())} Beatriz/literatura).")
    for fn in (analyze_seqfold, analyze_sizes, analyze_diversity,
               analyze_rarefaction, compute_descriptors):
        try:
            fn(df)
        except Exception as e:
            md(f"\n_⚠ {fn.__name__} falhou: {e}_")
    (ANADIR / "analysis_summary.md").write_text("\n".join(_md_lines), encoding="utf-8")
    print(f"\n✔ Resumo: {ANADIR / 'analysis_summary.md'}")
    print(f"✔ Figuras: {FIGDIR}")

if __name__ == "__main__":
    main()
