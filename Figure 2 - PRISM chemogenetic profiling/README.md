# PRISM chemogenetic profiling of XPR1 and IP6K inhibitors

Code and data for the PRISM panels of Figure 2 (2C-F) and Extended Data Figure 3 (3B-F) of
*Dynamic inositol pyrophosphate synthesis is a targetable therapeutic opportunity in ovarian
cancer* (Bondeson et al.).

475 barcoded cancer cell lines were pooled and treated with the XPR1 inhibitor RBD27, the IP6K
inhibitor SC919 or vehicle, and the barcodes were counted by sequencing (PRISM, Broad Institute).
The counts were processed with the CMAP **Sushi** pipeline
(<https://github.com/theprismlab/sushi>, commit `5be4ec2`, default settings). The Figshare
dataset holds the inputs and every intermediate table of that pipeline, from the raw barcode
counts to the log2 fold changes, so the analysis here starts at Sushi's output:
`figure2_analysis.Rmd` computes the sensitivity (AUC) of each cell line, relates it to
mRNA expression and gene dependency from DepMap Public 25Q3, regenerates every panel, writes a
source-data workbook, and prints each recomputed number next to its published value.

| Output file in `results/` | Panel |
|---|---|
| `fig2d_RBD27_vs_SC919_AUC.pdf` | Fig 2D: RBD27 against SC919 sensitivity, by SLC34A2 expression and lineage |
| `fig2e_SLC34A2_stratified_violin.pdf` | Fig 2E: sensitivity in each group, Mann-Whitney U tests |
| `fig2f_biomarker_concordance.pdf` | Fig 2F: correlation of expression and gene dependency with sensitivity, genome-wide |
| `ed3b_pathway_gene_Chronos_distributions.pdf` | ED Fig 3B: Chronos scores of the pathway genes across DepMap |
| `ed3c_cell_line_recovery.pdf` | ED Fig 3C: cell lines missing from each PCR well |
| `ed3d_expression_vs_sensitivity_logrank.pdf` | ED Fig 3D: expression and fold change across cell lines ranked by sensitivity |
| `ed3e_dependency_vs_sensitivity.pdf` | ED Fig 3E: correlations among the pathway genes and with sensitivity |
| `ed3f_PRISM_drug_drug_correlation.pdf` | ED Fig 3F: correlation with other PRISM compounds |
| `PRISM_sensitivity_by_cell_line.csv` | Fold change at each dose and AUC of the 475 cell lines |
| `Figure2_source_data.xlsx` | The data behind each panel and the statistics (sheet guide below) |
| `checks_vs_paper.csv` | Every recomputed number next to its published value |

> **Note on ED Fig 3D and 3F in this submission.** Two things will be corrected on re-submission, and
> the panels in `results/` already show the corrected versions.
> (1) The XPR1 mRNA correlations printed beside ED Fig 3D (-0.05 for RBD27, -0.057 for SC919) are
> replaced by the recomputed values, -0.052 and -0.059. The other 24 printed numbers agree to three decimals.
> (2) The title of ED Fig 3F changes. The panel shows how weakly the other compounds of the PRISM
> repurposing screen correlate with RBD27 and SC919 (all correlations below 0.3, the largest 0.19), against
> the correlation of 0.81 between RBD27 and SC919 themselves (Fig 2D). The submitted title, "RBD27:SC919
> correlation = 0.81", does not say that this is a different correlation from the ones plotted. The axes,
> labelled "Spearman" in the submitted panel, show Pearson correlations.

Not generated here: Fig 2A and 2C (schematics), Fig 2B (viability dose responses, with Fig 4F-G),
ED Fig 2 (protein engineering) and ED Fig 3A (flow cytometry).

## Get the data

The analysis reads the Figshare dataset *2026 IP6K paper - Figure 2 PRISM screen*
(<https://doi.org/10.6084/m9.figshare.34252608>; 66 MB). Download it into `data/` (see `data/README.md`), or keep it elsewhere and point to it:

```bash
export FIG2_DATA_DIR=/path/to/figshare/folder
```

## Run it

```bash
Rscript -e 'rmarkdown::render("figure2_analysis.Rmd")'
```

This needs pandoc (bundled with RStudio). Without pandoc, run the code directly:

```bash
Rscript -e 'knitr::knit("figure2_analysis.Rmd")'
```

Runtime is about 20 seconds. Every figure and the workbook are overwritten in `results/`. The
committed copies in `results/` are the expected output.

Tested with R 4.4.1 on macOS (arm64) and tidyverse 2.0.0 (dplyr 1.2.0, tidyr 1.3.1, readr 2.1.6,
ggplot2 4.0.2), here 1.0.2, scales 1.4.0, ggrepel 0.9.6, ggrastr 1.0.2, ggdist 3.3.3, patchwork 1.3.2,
broom 1.0.7, writexl 1.5.4, ragg 1.5.0, knitr 1.49. Install with:

```r
install.packages(c("tidyverse", "here", "scales", "ggrepel", "ggrastr", "ggdist", "patchwork",
                   "broom", "writexl", "ragg", "knitr", "rmarkdown"))
```

`.here` marks the project root for `here::here()`. Run from this folder and do not delete it.

## Layout

```
figure2_analysis.Rmd   the whole analysis
data/                  put the Figshare files here (README.md lists them)
results/               figures, source-data workbook, checks, expected output
```

## Analysis

### From Sushi output to sensitivity

1. **Check of the Sushi step.** The fold change of each cell line is recomputed from
   `normalized_counts.csv` (average of the three technical replicate wells, divided by the median
   of the vehicle wells of the same pool, log2) and compared with `l2fc.csv`: all 3,339 values agree
   to 1e-14. The normalization with the spiked-in control barcodes is Sushi's own step and is not
   repeated.
2. **Fold change and AUC.** The fold change of a cell line at a dose is Sushi's `median_l2fc`
   (`collapsed_l2fc.csv`; there is one biological replicate). The **AUC** of a drug is the mean over its
   doses of `min(2^l2fc, 1)`: it is 1 for a cell line that is not affected and approaches 0 for one that
   is eliminated at every dose. RBD27 has three doses (220, 670, 2000 nM) and SC919 four (100, 300, 1000,
   3000 nM). The AUC table is `results/PRISM_sensitivity_by_cell_line.csv`; it reproduces the AUC table the
   figures were first drawn from for all 475 cell lines (largest difference 5e-10).

### Panels

- **ED Fig 3C.** From `filtered_counts.csv`: for every PCR well, the number of cell lines with no reads
  (reported as missing) in each of the two pools. At most 5 of the 236 and 241 cell lines are missing from a well.
- **Fig 2D and 2E.** The 472 cell lines with an expression profile are split by lineage (ovary/fallopian tube
  and uterus together, against all others) and by SLC34A2 mRNA (high above 4 log2(TPM+1)). 2D plots RBD27
  against SC919 AUC, with the Pearson correlation; 2E shows each group's AUC and the two-sided Mann-Whitney U
  test of each group against the reference group (other lineages, SLC34A2-low).
- **Fig 2F.** The Pearson correlation of every gene's mRNA expression, and of every gene's Chronos score, with
  the AUC of each drug across the PRISM cell lines (genes measured in at least 100 of them: 19,203 expression
  and 17,931 dependency profiles). The PRISM biomarker module returns the same correlations for the
  top-ranked genes only; all 9,496 of those correlations are reproduced (largest difference 6e-6).
- **ED Fig 3B.** Chronos scores of 17 pathway genes in all 1,186 DepMap models (NUDT4: 45).
- **ED Fig 3D.** The cell lines, ranked by combined sensitivity, on a log-rank axis, with the mRNA expression
  of six genes, the AUC and the fold change at each dose. The numbers at the right are the Pearson
  correlations of each row with RBD27 (orange) and SC919 (teal) AUC.
- **ED Fig 3E.** Pearson correlations among the 17 pathway genes' Chronos scores and the two AUCs, in the 400
  PRISM cell lines with CRISPR data. Below, the median and variance of each profile with the sign chosen so that
  larger always means more killing (AUC is shown as 1 - AUC, Chronos is negated).
- **ED Fig 3F.** Correlations of RBD27 and SC919 AUC with the AUC of the 284 compounds of the public PRISM
  repurposing screen that were measured in at least 430 of the cell lines, taken from the PRISM biomarker output
  (the repurposing data themselves are not in the dataset). The five compounds with the highest and the five with
  the lowest mean correlation are labelled.

## Source data workbook

`results/Figure2_source_data.xlsx`. The `README` sheet inside lists each sheet; in brief:

| Sheet | Contents |
|---|---|
| `Sensitivity_by_cell_line` | 475 cell lines: lineage, fold change at each dose, AUC of each drug |
| `Fig2D_E_cell_lines`, `Fig2E_stats` | The 472 cell lines with SLC34A2 expression and group; group sizes, medians, U-test p-values |
| `Fig2F_correlations`, `Fig2F_vs_biomarker_output` | Correlation of every gene with each drug; comparison with the PRISM biomarker output |
| `ED3B_Chronos`, `ED3C_recovery`, `ED3D_correlations`, `ED3E_correlations`, `ED3F_repurposing` | The data behind each Extended Data panel |
| `Sushi_l2fc_check` | Fold changes recomputed from the normalized counts next to Sushi's |
| `checks_vs_paper` | Every recomputed number next to its published value |

## Checks

`results/checks_vs_paper.csv` lists 54 numbers; 52 reproduce. They include the 475 cell lines, the 472 with
expression, the group sizes of Fig 2D (16, 12, 36, 408), R = 0.81 and the six U-test p-values (1.1e-5, 0.520,
0.003, 0.003, 0.721, 0.016) of Fig 2E, the pool sizes and the maximum of 5 missing cell lines of ED Fig 3C,
the 26 correlations printed beside ED Fig 3D, the 1,186 and 45 cell lines of ED Fig 3B, and the ten
compounds labelled in ED Fig 3F, whose correlations with either drug are all below 0.3 (largest 0.19).
The two that do not match are the XPR1 mRNA correlations of ED Fig 3D (see the note at the top).
