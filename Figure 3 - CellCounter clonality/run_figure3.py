#!/usr/bin/env python3
"""Figure 3D-F and Extended Data Figure 4G: CellCounter clonality of OVISE xenografts.

Reads the Figshare dataset (see data/README.md), regenerates the CellCounter panels,
writes a source-data workbook, and prints a table that compares each number quoted in the
paper with the value recomputed here.

    python run_figure3.py --data-dir data --out-dir results

Panels: Fig 3D (ascites bioluminescence against CellCounter reads), 3E (clone size by rank),
3F (unique clones against Shannon evenness) and ED Fig 4G (clone overlap between samples).
"""
import argparse
import os
import string
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import stats  # noqa: E402

FILES = {
    "counts": "OVISE_xenograft_CellCounter_UMI_read_counts.csv.gz",
    "animals": "OVISE_IP_xenograft_animals_randomization_and_outcome.csv",
    "ex_vivo": "OVISE_IP_xenograft_ascites_ex_vivo_BLI.csv",
}

TEAL, GREY = "darkcyan", "#808080"
SC919, VEH, PRE = "SC919 25 mg/kg QD", "Vehicle", "Pre-injection"
GROUP_COLOR = {SC919: TEAL, VEH: GREY, PRE: "sienna"}
ORGANS = ["ascites", "kidney", "liver", "lung", "ovary"]
ORGAN_COLOR = {"ascites": "#E69F00", "kidney": "#56B4E9", "liver": "#009E73",
               "lung": "#F0E442", "ovary": "#CC79A7", "preinj": "#999999"}
CHECKS = []   # (panel, quantity, recomputed, expected, expected_source, ok)


# ---------------------------------------------------------------------------------------
# setup
# ---------------------------------------------------------------------------------------
def set_style():
    mpl.rcParams["pdf.fonttype"] = 42      # editable text in Illustrator
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    mpl.rcParams["axes.linewidth"] = 1.5
    mpl.rcParams["xtick.major.width"] = 1.5
    mpl.rcParams["ytick.major.width"] = 1.5
    mpl.rcParams["font.size"] = 9


def clean_axes(ax):
    sns.despine(ax=ax, top=True, right=True)


def save(fig, out, name):
    fig.savefig(Path(out) / name, bbox_inches="tight", transparent=True)
    plt.close(fig)


def check(panel, quantity, got, expected, source, tol=None, text=False):
    if text:
        ok = str(got) == str(expected)
    elif tol is None:
        ok = bool(np.isclose(got, expected, rtol=1e-3))
    else:
        ok = bool(abs(got - expected) <= tol)
    CHECKS.append((panel, quantity, got, expected, source, ok))


# ---------------------------------------------------------------------------------------
# data loading
# ---------------------------------------------------------------------------------------
def load_counts(path, cache=None):
    """UMI read counts of the 75 samples shown in the paper: 12 IP mice x 5 organs and the 15
    pre-injection samples (the Figshare table holds only these)."""
    if cache and Path(cache).exists():
        return pd.read_pickle(cache)
    use = ["bio_rep", "organ", "umi", "n", "n_umis", "umi_rank"]
    parts = []
    for ch in pd.read_csv(path, usecols=use, chunksize=1_000_000, dtype={"bio_rep": str, "organ": str, "umi": str}):
        ch["organ_b"] = ch.organ.str.split("-").str[-1].str.strip().str.lower()
        parts.append(ch)
    df = pd.concat(parts, ignore_index=True)
    df["sample"] = df.bio_rep + "_" + df.organ
    for c in ["bio_rep", "organ", "organ_b", "sample"]:
        df[c] = df[c].astype("category")
    assert df["sample"].nunique() == 75, df["sample"].nunique()
    if cache:
        df.to_pickle(cache)
    return df


def load_tables(data_dir):
    d = Path(data_dir)
    return {k: pd.read_csv(d / v) for k, v in FILES.items() if k != "counts"}


# ---------------------------------------------------------------------------------------
# CellCounter metrics
# ---------------------------------------------------------------------------------------
def sample_table(counts, animals):
    """One row per sequenced sample with group, reads, clones and Shannon evenness."""
    g = counts.groupby("sample", observed=True)
    s = g.agg(bio_rep=("bio_rep", "first"), organ=("organ", "first"), organ_b=("organ_b", "first"),
              n_umis=("n_umis", "first"), clones_in_file=("n", "size"), total_reads=("n", "sum")).reset_index()
    # Shannon diversity H' = -sum p ln p over the clones in the file; evenness E = H'/ln(S), S = clones in the file
    cn = counts[["sample", "n"]].copy()
    cn["p"] = cn.n / cn.groupby("sample", observed=True).n.transform("sum")
    cn["plogp"] = cn.p * np.log(cn.p)
    H = -cn.groupby("sample", observed=True).plogp.sum()
    S = cn.groupby("sample", observed=True).size()
    s["shannon_index"] = s["sample"].map(H)
    s["shannon_evenness"] = s["sample"].map(H / np.log(S))
    gmap = animals.dropna(subset=["mouse_id_pipeline"]).set_index("mouse_id_pipeline").group
    s["group"] = np.where(s.organ_b == "preinj", PRE, s.bio_rep.astype(str).map(gmap))
    s["mouse_id_pipeline"] = np.where(s.organ_b == "preinj", None, s.bio_rep.astype(str))
    return s


def overlap_matrix(counts, samples, top_fraction=0.1):
    """Overlap coefficient |A and B| / min(|A|, |B|) between the top 10% of clones (by reads)
    of every pair of samples, as used for ED Fig 4G."""
    sub = counts[counts["sample"].isin(samples)]
    sets = {}
    for s, g in sub.groupby("sample", observed=True):
        k = int(len(g) * top_fraction)
        sets[s] = set(g.loc[g.umi_rank < k, "umi"])
    names = sorted(sets)
    m = pd.DataFrame(np.nan, index=names, columns=names)
    for a in range(len(names)):
        for b in range(a, len(names)):
            A, B = sets[names[a]], sets[names[b]]
            v = len(A & B) / min(len(A), len(B)) if A and B else 0.0
            m.iat[a, b] = m.iat[b, a] = v
    return m


# ---------------------------------------------------------------------------------------
# Fig 3D, 3E, 3F: CellCounter
# ---------------------------------------------------------------------------------------
def fig3d(samples, ev, out):
    a = samples[(samples.organ_b == "ascites") & samples.group.isin([VEH, SC919])].merge(
        ev[["mouse_id_pipeline", "avg_radiance"]], on="mouse_id_pipeline")
    assert len(a) == 12, len(a)
    r, p = stats.spearmanr(a.total_reads, a.avg_radiance)
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    for g in [VEH, SC919]:
        d = a[a.group == g]
        ax.scatter(np.log10(d.total_reads), np.log10(d.avg_radiance), s=45, color=GROUP_COLOR[g], edgecolor="black", linewidth=0.9, zorder=3)
    ax.text(0.96, 0.05, f"Spearman R: {r:.2f}\np-value: {p:.2e}", transform=ax.transAxes, ha="right", va="bottom",
            bbox=dict(boxstyle="square,pad=0.4", fc="white", ec="black"))
    ax.set_xlim(4, 6.7)
    ax.set_ylim(7, 9.2)
    ax.set_xlabel("CellCounter Total Reads, Log10", fontweight="bold")
    ax.set_ylabel("Bioluminescence, Log10", fontweight="bold")
    ax.set_title("Ascites tumor burden", fontweight="bold", fontsize=9)
    clean_axes(ax)
    save(fig, out, "fig3d_ascites_BLI_vs_CellCounter_reads.pdf")
    check("Fig 3D", "Spearman R", round(r, 2), 0.72, "Fig 3D panel and legend")
    check("Fig 3D", "Spearman p", round(p, 5), 0.00824, "Fig 3D panel (the legend prints 6.24e-3)", tol=5e-6)
    return a.assign(log10_total_reads=np.log10(a.total_reads), log10_bli=np.log10(a.avg_radiance)), r, p


def fig3e(counts, samples, out):
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    sel = samples[samples.organ_b.isin(["ascites", "preinj"]) & samples.group.isin([VEH, SC919, PRE])]
    pal = {VEH: sns.color_palette("Greys", n_colors=(sel.group == VEH).sum() + 4)[4:],
           SC919: sns.light_palette("teal", n_colors=(sel.group == SC919).sum() + 3)[3:],
           PRE: sns.light_palette("sienna", n_colors=(sel.group == PRE).sum() + 3)[3:]}
    by_sample = {s: g for s, g in counts[counts["sample"].isin(sel["sample"])].groupby("sample", observed=True)}
    clone_rows = []
    for z, grp in enumerate([PRE, VEH, SC919], start=1):
        names = sorted(sel[sel.group == grp]["sample"])
        for i, s in enumerate(names):
            g = by_sample[s]
            x, y = np.log10(g.umi_rank.values + 1), np.log10(g.n.values)
            ax.scatter(x, y, color=[pal[grp][i]], s=6, alpha=0.7, edgecolors="none", zorder=z, rasterized=True)
            if grp != PRE:
                clone_rows.append(pd.DataFrame({"sample": s, "group": grp, "clone_rank": g.umi_rank.values + 1, "reads": g.n.values}))
    ax.set_xlim(-0.3, 6)
    ax.set_ylim(0.3, 5.3)
    ax.set_xlabel("CellCounter Clone rank, log10", fontweight="bold")
    ax.set_ylabel("CellCounter clone size\n(number of reads, Log10)", fontweight="bold")
    clean_axes(ax)
    save(fig, out, "fig3e_clone_size_by_rank.pdf")
    return pd.concat(clone_rows, ignore_index=True)


def fig3f(samples, out):
    s = samples[samples.organ_b.isin(["ascites", "preinj"]) & samples.group.isin([VEH, SC919, PRE])].copy()
    s["log10_unique_clones"] = np.log10(s.n_umis)
    fig, ax = plt.subplots(figsize=(2.4, 2.2))
    for g in [PRE, VEH, SC919]:
        d = s[s.group == g]
        ax.scatter(d.log10_unique_clones, d.shannon_evenness, s=38, color=GROUP_COLOR[g], edgecolor="black", linewidth=0.8, alpha=0.9, zorder=3)
    ax.text(0.3, 0.88, "SC919", transform=ax.transAxes, color=TEAL, fontweight="bold", fontsize=8)
    ax.text(0.6, 0.55, "Vehicle", transform=ax.transAxes, color=GREY, fontweight="bold", fontsize=8)
    ax.text(0.55, 1.02, "Pre-injection", transform=ax.transAxes, color="sienna", fontweight="bold", fontsize=8)
    ax.set_xlim(3, 6)
    ax.set_ylim(0.55, 1.03)
    ax.set_xlabel("Unique clones (log10)", fontweight="bold")
    ax.set_ylabel("Shannon Evenness", fontweight="bold")
    clean_axes(ax)
    save(fig, out, "fig3f_unique_clones_vs_evenness.pdf")
    med = s[s.group.isin([VEH, SC919])].groupby("group").n_umis.median()
    check("Fig 3F / Results", "median unique clones, vehicle", med[VEH], 26000, "Results text (~26,000)", tol=500)
    check("Fig 3F / Results", "median unique clones, SC919", med[SC919], 2300, "Results text (~2,300)", tol=100)
    check("Fig 3F / Results", "fold difference in unique clones", round(med[VEH] / med[SC919], 1), 11, "Results text (~11-fold)", tol=0.5)
    return s


# ---------------------------------------------------------------------------------------
# ED4G: clone overlap between samples
# ---------------------------------------------------------------------------------------
def ed4g(counts, samples, out):
    treated = samples[samples.organ_b.isin(ORGANS) & samples.group.isin([VEH, SC919])]["sample"].tolist()
    pre = samples[samples.organ_b == "preinj"]["sample"].tolist()
    mice = sorted(samples[samples.group.isin([VEH, SC919])].bio_rep.astype(str).unique())
    sc = [m for m in mice if samples[(samples.bio_rep.astype(str) == m)].group.iloc[0] == SC919]
    ve = [m for m in mice if m not in sc]
    label = {m: f"SC{i + 1}" for i, m in enumerate(sc)} | {m: f"V{i + 1}" for i, m in enumerate(ve)}
    label |= {f"M24-X11_PREINJ-D4-{c}": c for c in string.ascii_uppercase}
    mt = overlap_matrix(counts, treated)
    mp = overlap_matrix(counts, pre)
    for m, trig, name in [(mt, "lower", "ed4g_overlap_treated_lower.pdf"), (mp, "upper", "ed4g_overlap_preinjection_upper.pdf")]:
        draw_overlap(m, trig, label, out, name)
    return mt, mp, label


def draw_overlap(m, triangle, label, out, name):
    n = len(m)
    organs = [next((o for o in ORGAN_COLOR if o in s.lower()), "unknown") for s in m.index]
    bio = [label.get(s.rsplit("_", 1)[0] if "PREINJ" in s else s.split("_IP")[0], s) for s in m.index]
    fig, ax = plt.subplots(figsize=(14, 12))
    mask = np.tril(np.ones((n, n), bool)) if triangle == "upper" else np.triu(np.ones((n, n), bool))
    sns.heatmap(m, mask=mask, cmap="viridis", square=True, ax=ax, vmin=0, vmax=1, cbar_kws={"label": "Overlap Coefficient", "shrink": 0.8})
    ax.set_xticks([])
    ax.set_yticks([])
    tw = 1.0
    if triangle == "upper":
        ax.set_xlim(0, n + tw)
        ax.set_ylim(n, -tw)
    else:
        ax.set_xlim(-tw, n)
        ax.set_ylim(n + tw, 0)
    for i, o in enumerate(organs):
        c = ORGAN_COLOR.get(o, "#FF00FF")
        if triangle == "upper":
            ax.bar(i + 0.5, tw, bottom=-tw, width=1.0, color=c, edgecolor="white", linewidth=0.5, zorder=5)
            ax.barh(i + 0.5, tw, left=n, height=1.0, color=c, edgecolor="white", linewidth=0.5, zorder=5)
        else:
            ax.bar(i + 0.5, tw, bottom=n, width=1.0, color=c, edgecolor="white", linewidth=0.5, zorder=5)
            ax.barh(i + 0.5, tw, left=-tw, height=1.0, color=c, edgecolor="white", linewidth=0.5, zorder=5)
    start = 0
    blocks = []
    for i in range(1, n):
        if bio[i] != bio[i - 1]:
            blocks.append((bio[i - 1], start, i))
            start = i
    blocks.append((bio[-1], start, n))
    for rep, a, b in blocks:
        mid = (a + b) / 2
        if a > 0:
            lw, col = (3, "black") if bio[a - 1][:2] != rep[:2] else (1, "white")
            if triangle == "upper":
                ax.hlines(a, a, n + tw, color=col, linewidth=lw, zorder=10)
                ax.vlines(a, -tw, a, color=col, linewidth=lw, zorder=10)
            else:
                ax.hlines(a, -tw, a, color=col, linewidth=lw, zorder=10)
                ax.vlines(a, a, n + tw, color=col, linewidth=lw, zorder=10)
        if triangle == "upper":
            ax.text(n + tw + 0.5, mid, rep, rotation=-90, va="center", ha="left", fontsize=14, fontweight="bold")
            ax.text(mid, -tw - 0.5, rep, ha="center", va="bottom", fontsize=14, fontweight="bold")
        else:
            ax.text(-tw - 0.5, mid, rep, rotation=90, va="center", ha="center", fontsize=14, fontweight="bold")
            ax.text(mid, n + tw + 0.5, rep, ha="center", va="top", fontsize=14, fontweight="bold")
    handles = [Patch(facecolor=c, edgecolor="black", label=o.capitalize()) for o, c in ORGAN_COLOR.items() if o in set(organs)]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.35 if triangle == "upper" else 1.25, 1), title="Organ Key", frameon=False, fontsize=12, title_fontsize=14)
    fig.savefig(Path(out) / name, bbox_inches="tight", transparent=True)
    plt.close(fig)


# ---------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=os.environ.get("FIG3_DATA_DIR", "data"), help="folder with the Figshare files")
    ap.add_argument("--out-dir", default="results")
    ap.add_argument("--cache", default=None, help="optional pickle that caches the count table between runs")
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    set_style()
    t0 = time.time()

    t = load_tables(a.data_dir)
    counts = load_counts(Path(a.data_dir) / FILES["counts"], a.cache)
    print(f"counts loaded: {len(counts):,} clone rows, {counts['sample'].nunique()} samples ({time.time() - t0:.0f}s)")
    animals = t["animals"]
    samples = sample_table(counts, animals)

    # animal accounting
    ng = animals.groupby("group").size()
    check("Animals", "IP mice randomized (SC919 / vehicle)", f"{ng[SC919]} / {ng[VEH]}", "8 / 7", "Randomization_IP_d14 sheet", text=True)
    ok = animals[animals.outcome == "analysed at endpoint"].groupby("group").size()
    check("Animals", "IP mice analysed (SC919 / vehicle)", f"{ok[SC919]} / {ok[VEH]}", "5 / 7", "Fig 3C legend", text=True)

    d_tab, rho, pval = fig3d(samples, t["ex_vivo"], out)
    clones = fig3e(counts, samples, out)
    f_tab = fig3f(samples, out)
    print(f"figures drawn ({time.time() - t0:.0f}s)")
    mt, mp, label = ed4g(counts, samples, out)
    print(f"overlap heatmaps drawn ({time.time() - t0:.0f}s)")
    check("ED Fig 4G", "samples in treated heatmap (12 mice x 5 organs)", len(mt), 60, "ED Fig 4G")
    check("ED Fig 4G", "pre-injection samples", len(mp), 15, "ED Fig 4G")

    # sensitivity of the clone counts to the read threshold used for 'unique clones'
    sens = []
    for thr in (3, 10):
        k = counts[(counts.n >= thr) & counts["sample"].isin(f_tab[f_tab.group.isin([VEH, SC919])]["sample"])].groupby("sample", observed=True).size()
        m = f_tab[f_tab.group.isin([VEH, SC919])].set_index("sample").assign(k=k).groupby("group").k.median()
        sens.append((f"clones with >= {thr} reads, counted from the count table", m[VEH], m[SC919], m[VEH] / m[SC919]))
    m = f_tab[f_tab.group.isin([VEH, SC919])].groupby("group").n_umis.median()
    sens.insert(0, ("n_umis column (threshold not recorded; used for Fig 3F and the Results text)", m[VEH], m[SC919], m[VEH] / m[SC919]))
    sens = pd.DataFrame(sens, columns=["definition_of_unique_clones", "median_vehicle", "median_SC919", "fold_difference"])

    chk = pd.DataFrame(CHECKS, columns=["panel", "quantity", "recomputed", "expected", "expected_source", "match"])
    pd.set_option("display.width", 220, "display.max_colwidth", 70, "display.max_rows", 200)
    print("\nRecomputed values against the paper\n")
    print(chk.to_string(index=False))
    print("\nUnique-clone fold difference by definition (Methods state clones with < 10 reads were excluded)\n")
    print(sens.round(1).to_string(index=False))

    guide = pd.DataFrame([
        ("Fig3D", "Ascites BLI against CellCounter total reads, 12 mice"),
        ("Fig3E_ascites_clones", "Every clone of the 12 ascites samples (rank and reads). Pre-injection clones (4.75 million rows) are in the Figshare counts file"),
        ("Fig3F", "Unique clones and Shannon evenness for 27 samples"),
        ("EDFig4G_*", "Overlap coefficients (top 10% of clones) between samples"),
        ("unique_clone_sensitivity", "Median unique clones under three definitions"),
        ("checks_vs_paper", "Recomputed value against the paper for every quoted number"),
    ], columns=["sheet", "contents"])
    with pd.ExcelWriter(out / "Figure3_CellCounter_source_data.xlsx", engine="openpyxl") as xw:
        guide.to_excel(xw, sheet_name="README", index=False)
        d_tab[["mouse_id_pipeline", "group", "total_reads", "avg_radiance", "log10_total_reads", "log10_bli"]].rename(
            columns={"avg_radiance": "ascites_avg_radiance"}).to_excel(xw, sheet_name="Fig3D", index=False)
        clones.to_excel(xw, sheet_name="Fig3E_ascites_clones", index=False)
        f_tab[["sample", "group", "organ_b", "n_umis", "clones_in_file", "total_reads", "shannon_index", "shannon_evenness"]].to_excel(xw, sheet_name="Fig3F", index=False)
        mt.to_excel(xw, sheet_name="EDFig4G_treated")
        mp.to_excel(xw, sheet_name="EDFig4G_preinjection")
        sens.to_excel(xw, sheet_name="unique_clone_sensitivity", index=False)
        chk.to_excel(xw, sheet_name="checks_vs_paper", index=False)
    chk.to_csv(out / "checks_vs_paper.csv", index=False)
    n_bad = int((~chk.match).sum())
    print(f"\n{len(chk) - n_bad} of {len(chk)} checks match.")
    print(f"done in {time.time() - t0:.0f}s; outputs in {out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
