# Genome-scale modifier screens of XPR1 and IP6K inhibition

Code and data for Figure 4A and Extended Data Figure 5B of *Dynamic inositol pyrophosphate
synthesis is a targetable therapeutic opportunity in ovarian cancer* (Bondeson et al.).

The Humagne enCas12a library (two sets, C and D, of about 20,000 constructs each) was screened in
RMGI ovarian cancer cells for 14 days with vehicle, 200 nM RBD27 or 350 nM SC919, and the
constructs were counted with PoolQ. Starting from the summed read counts of the two sets,
`figure4_analysis.Rmd` computes the fold change of every construct, scores every gene with the
hypergeometric method of the Broad GPP portal, regenerates the two panels, writes a source-data
workbook, and checks the gene scores against the portal's own output.

| Output file in `results/` | Panel |
|---|---|
| `fig4a_modifier_screen_RBD27_vs_SC919.pdf` | Fig 4A: gene-level rescue by RBD27 and SC919, coloured by RMGI Chronos score |
| `ed5b_screen_sample_concordance.pdf` | ED Fig 5B: agreement between the screen samples and with RMGI Chronos |
| `Figure4_modifier_screen_source_data.xlsx` | Counts, construct and gene scores and the data behind each panel (sheet guide below) |
| `checks_vs_paper.csv` | Every recomputed number next to its published or portal value |

The other panels of Figure 4 and Extended Data Figure 5 (competition assay, phosphate efflux,
PiT1 and viability dose responses, immunoblot and histology, screen growth curves, RMGI xenografts)
are not generated here.

## Get the data

The analysis reads the Figshare dataset *2026 IP6K paper - Figure 4 genome-scale drug modifier screen*
(<https://doi.org/10.6084/m9.figshare.34253274>; 7 data files, 21 MB). Download it into `data/` (see `data/README.md` for the
file list), or keep it elsewhere and point to it:

```bash
export FIG4_DATA_DIR=/path/to/figshare/folder
```

## Run it

```bash
Rscript -e 'rmarkdown::render("figure4_analysis.Rmd")'
```

This needs pandoc (bundled with RStudio). Without pandoc, run the code directly:

```bash
Rscript -e 'knitr::knit("figure4_analysis.Rmd")'
```

Runtime is about 20 seconds. Every figure and the workbook are overwritten in `results/`. The
committed copies in `results/` are the expected output.

Tested with R 4.4.1 on macOS (arm64) and tidyverse 2.0.0 (dplyr 1.2.0, tidyr 1.3.1, readr 2.1.6,
ggplot2 4.0.2), scales 1.4.0, ggrepel 0.9.6, scattermore 1.2, reshape2 1.4.5, janitor 2.2.1,
writexl 1.5.4, here 1.0.2, knitr 1.49. Install with:

```r
install.packages(c("tidyverse", "scales", "ggrepel", "scattermore", "reshape2", "janitor",
                   "writexl", "here", "knitr", "rmarkdown"))
```

`.here` marks the project root for `here::here()`. Run from this folder and do not delete it.

## Layout

```
figure4_analysis.Rmd   the whole analysis
data/                  put the Figshare files here (README.md lists them)
results/               figures, source-data workbook, checks, expected output
```

## Analysis

Each library set has its own count table (`RMGI_Humagne_Set{C,D}_summed_read_counts.csv`).

1. **Counts to log2 RPM.** `log2(count / sample total x 10^6 + 1)`, within each sample of a set.
2. **Construct fold change.** Within each set, the median log2 RPM of the drug replicates
   minus the median log2 RPM of the DMSO replicates. Set D has one RBD27 replicate.
3. **Gene scores (Broad GPP portal method, reimplemented).** All 40,376 constructs of both
   sets are ranked together by fold change. For a gene with n constructs, the k-th best
   construct, at rank r, gets p = P(exactly k of the gene's n constructs fall among the r
   lowest-ranked constructs), the hypergeometric probability with N = 40,376, once from
   each end of the ranking. The gene's score is the larger of the two mean -log10(p), and
   its average LFC is the mean of its constructs' LFC. Fig 4A plots the average LFC.
4. **Check against the portal.** The portal's own output files are in the dataset. For the
   19,752 genes with the same constructs in both, the recomputed average LFC and score
   agree with the portal to 1e-10 for both drugs. The 599 control pseudo-genes
   (`NO_SITE_*`, `ONE_INTERGENIC_SITE_*`) pair control constructs two at a time in the
   portal; the pairing is not recoverable from the design file, so they are not reproduced
   and are not plotted.
5. **Fig 4A.** Average LFC for RBD27 (x) against SC919 (y), coloured by RMGI Chronos
   (DepMap 25Q3); the genes labelled in the paper are highlighted. The three correlations
   printed on the panel (-0.100, -0.037, 0.295) are recomputed.
6. **ED Fig 5B.** Gene-level log2 fold change of each sample against `Presort` (constructs
   with fewer than 5 early reads, and control constructs, left out) for the 11 samples of
   the two sets, with RMGI Chronos; Pearson correlations of the lower triangle.

The Methods say that arrays with fewer than 5 reads in the early sample were excluded from
Fig 4A. The portal run did not exclude any: every construct fold change regenerates from the
unfiltered counts, so Fig 4A is computed on all constructs, as the published values were.

## Source data workbook

`results/Figure4_modifier_screen_source_data.xlsx`. The `README` sheet inside lists each sheet:

| Sheet | Contents |
|---|---|
| `Counts_SetC`, `Counts_SetD`, `Sample_key` | The input count tables and how each column was assembled |
| `Fig4A_genes` | Gene-level average LFC and score for both drugs, with RMGI Chronos |
| `Fig4A_constructs` | Construct-level log2 RPM and LFC against DMSO for all 40,376 constructs |
| `Fig4A_portal_check` | Recomputed against portal output, per drug |
| `ED5B_data`, `ED5B_correlation` | Gene-level fold changes per sample; the correlation matrix |
| `checks_vs_paper` | Recomputed numbers next to published ones |

## Checks

`results/checks_vs_paper.csv` lists 7 numbers and all reproduce: the gene scores of both drugs against
the portal output (largest difference 7e-11) and the three correlations printed on Fig 4A.
