# Base editor screen of XPR1 and KIDINS220 in OVISE cells

Code and data for the base editor screen in Figure 1 and Extended Data Figure 1 of
*Dynamic inositol pyrophosphate synthesis is a targetable therapeutic opportunity
in ovarian cancer* (Bondeson et al.).

Starting from the raw sgRNA read counts, `base_editor_screen_analysis.Rmd`
reproduces every screen-derived panel and writes a source-data workbook with the
raw counts and the data behind each panel.

> **Note on ED Fig 1E in this submission.** The scoring counts on ED Fig 1E in the
> submitted manuscript (A>G 176 / 323 / 299; C>T 183 / 257 / 99) were computed from
> the **mean** of the three replicate log2 fold changes. Every other panel, the 211
> scoring sgRNA and Supplementary Table 1 use the **median**. **Upon re-submission
> we will switch ED Fig 1E to the median**, which gives A>G 203 / 335 / 322 and
> C>T 181 / 257 / 110.
>
> `results/edfig1e_scoring_distribution.pdf` is the median version that will
> replace the submitted panel. The mean-based numbers in the submitted panel are
> reproduced in the `EDFig1E_mean_of_reps` sheet of the source-data workbook.

| Output file in `results/` | Panel |
|---|---|
| `fig1a_scores_by_position.pdf` | Fig 1A: sgRNA depletion along XPR1 and KIDINS220 |
| `fig1b_domain_enrichment.pdf` | Fig 1B: domain-level enrichment of depleted sgRNA |
| `edfig1a_coverage.pdf` | ED Fig 1A: residue coverage of the library |
| `edfig1b_library_complexity.pdf` | ED Fig 1B: unique edits vs. saturating mutagenesis |
| `edfig1c_replicate_correlation.pdf` | ED Fig 1C: sample-to-sample correlation |
| `edfig1d_lfc_vs_rpm.pdf` | ED Fig 1D: depletion vs. Day 4 representation |
| `edfig1e_scoring_distribution.pdf` | ED Fig 1E: scoring sgRNA per timepoint (median-based; differs from the submitted panel, see the note above) |
| `Figure1_base_editor_screen_source_data.xlsx` | Raw counts and per-panel data (sheet guide below) |
| `sgRNA_map_positions_and_domains.csv` | Amino acid position, domain and name assigned to every sgRNA |

Panels 1C-1F and ED Fig 1F (validation assays, structure modelling) are not
generated here.

## Get the data

The analysis reads the Figshare dataset *2026 IP6K paper - Figure 1 base editor screen*
(<https://doi.org/10.6084/m9.figshare.34252485>; 14 data files, 7 MB). Download it into `data/` (see `data/README.md`
for the file list), or keep it elsewhere and point to it:

```bash
export FIG1_DATA_DIR=/path/to/figshare/folder
```

## Run it

```bash
Rscript -e 'rmarkdown::render("base_editor_screen_analysis.Rmd")'
```

This needs pandoc (bundled with RStudio). Without pandoc, run the code directly:

```bash
Rscript -e 'knitr::knit("base_editor_screen_analysis.Rmd")'
```

Runtime is about 15 seconds. Every figure and the workbook are overwritten in
`results/`. The committed copies in `results/` are the expected output.

Tested with R 4.4.1 on macOS (arm64) and tidyverse 2.0.0 (dplyr 1.2.0, tidyr 1.3.1,
readr 2.1.6, ggplot2 4.0.2), reshape2 1.4.5, scales 1.4.0, ggrepel 0.9.6, writexl
1.5.4, here 1.0.2, knitr 1.49. Install with:

```r
install.packages(c("tidyverse", "reshape2", "scales", "ggrepel", "writexl",
                   "here", "knitr", "rmarkdown"))
```

`.here` marks the project root for `here::here()`. Run from this folder and do not
delete it.

## Layout

```
base_editor_screen_analysis.Rmd   the whole analysis
data/                             put the Figshare files here (README.md lists them)
results/                          figures, source-data workbook, expected output
```

The Figshare dataset holds:

```
base_editor_library_sgRNA_design_GPP.csv    library design: predicted edits for both editors
base_editor_library_gene_annotations.csv    gene -> XPR1 / KIDINS220 / Positive Control / Negative Control
poolq_outputs/                              PoolQ output for both plates (see below)
```

## Data

The library of 4,039 sgRNA (1,011 targeting XPR1, 2,435 KIDINS220, the rest
controls) was cloned into two vectors and screened in OVISE cells. The Broad
Genetic Perturbation Platform (GPP) sequenced the screens (batch 2902200,
December 2021) and counted reads with **PoolQ v3.3.4**. The exact PoolQ command
line is in each `runinfo-*` file.

| Plate | Vector | Editor | Label in the code | Samples |
|---|---|---|---|---|
| p2 (`NP_GPP3155`) | pRDA_479 | ABE8e, A>G | `A.G` | Rep A-C at Days 4, 7, 14, 21, plus pDNA |
| p1 (`NP_GPP3154`) | pRDA_478 | BE3.9, C>T | `C.T` | Rep A-C at Days 4, 7, 21; Rep A-B at Day 14; four pDNA aliquots |

Files in `poolq_outputs/`, one set per plate:

| File | Contents |
|---|---|
| `counts-*` | **Raw read counts**, one row per sgRNA, PCR replicates already summed. This is the input to the analysis. |
| `barcodecounts-*` | Counts per sample-index barcode, before samples are assembled with the conditions file. Not used by the analysis. |
| `lognorm-*` | PoolQ's own log2 reads-per-million. Used only to check the conversion from raw counts. |
| `quality-*`, `correlation-*`, `runinfo-*` | PoolQ run QC, replicate correlation and run settings. |

Column headers are PoolQ's. Columns named `please spike in pDNA...` are GPP's
placeholder names for the plasmid library. Plate 1 has four of them, sequenced at
different depths (1.5-11.8 million reads); they correlate at r >= 0.96 with each
other and are averaged where a Day 0 value is needed (ED Fig 1C only).

## Analysis

1. **Counts to log2 RPM.** `log2(count / sample total x 10^6 + 1)`. This matches
   PoolQ's `lognorm` file to within 1e-14, and the analysis stops if it does not.
2. **sgRNA mapping.** Each sgRNA is assigned the amino acid position of its
   predicted edit. A guide with no predicted edit takes the position of the nearest
   guide on the same strand (`aa_directly_predicted` records which is which).
   Domains: XPR1 SPX < 224, Core 224-440, EXS > 440 (called Channel in the paper
   figure); KIDINS220 Ank <= 500, TM 501-1200, Disordered > 1200.
3. **Depletion.** Log2 fold change of each replicate against the **mean Day 4**
   log2 RPM (the first post-selection timepoint, not the plasmid). An sgRNA's
   depletion is the **median** of its replicates. The scoring cut-off is the
   median minus 2 standard deviations of the negative-control sgRNA, computed per
   editor and per day.
4. **Domain enrichment (Fig 1B).** sgRNA targeting XPR1 or KIDINS220 are ranked by
   Day 21 A>G depletion. A GSEA-style running score, weighted by |LFC|, is computed
   for each domain. The normalised enrichment score is the maximum score divided
   by the mean maximum of 1,000 permutations (seed 42); p is the fraction of
   permutations at least as extreme. Disordered is depleted rather than enriched
   and every permutation maximum is positive, so its NES is undefined (NaN here,
   "Na" in the paper figure).

**Mean and median.** Fig 1A-B, the scoring sgRNA reported in the text and
Supplementary Table 1 use the median of replicates. The scoring counts printed on
the submitted ED Fig 1E were produced with the mean of replicates instead; ED Fig
1E will be switched to the median upon re-submission (see the note at the top).
The analysis computes both: `edfig1e_scoring_distribution.pdf` uses the median,
and the `EDFig1E_mean_of_reps` sheet of the workbook reproduces the submitted
panel.

## Source data workbook

`results/Figure1_base_editor_screen_source_data.xlsx`. The `README` sheet inside
lists each sheet; in brief:

| Sheet | Contents |
|---|---|
| `Counts_A>G_ABE8e`, `Counts_C>T_BE3.9` | Raw read counts (4,039 sgRNA each) |
| `sgRNA_annotation` | Gene, control class, predicted edit, mutation class, position, domain and name of every sgRNA, per editor |
| `Fig1A` | Day 21 depletion, cut-off and score for every plotted sgRNA |
| `Fig1B_curves`, `Fig1B_stats`, `Fig1B_cutoff` | Running enrichment curves, NES and p per domain, rank of the cut-off |
| `EDFig1A` ... `EDFig1E` | Data behind each Extended Data panel |
| `EDFig1E_mean_of_reps` | ED Fig 1E scored on replicate means |

## Checks

The last section of the analysis recomputes numbers quoted in the paper from the
raw counts and prints them next to the published values. Run from raw counts, it
reproduces the 4,039-sgRNA library, the 211 scoring A>G sgRNA at Day 21 (114 in
XPR1, 97 in KIDINS220), the ED Fig 1A and 1B totals, and all five NES values in
Fig 1B. The submitted ED Fig 1E counts are reproduced by the mean-of-replicates
calculation, not the median one. The 211 sgRNA, their depletion values and their
cut-offs are identical to the
submitted Supplementary Table 1.
