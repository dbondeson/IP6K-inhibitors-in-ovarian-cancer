# IP6K inhibitors in ovarian cancer: code for the figures

Code that regenerates the sequencing-based figure panels of *Dynamic inositol pyrophosphate synthesis
is a targetable therapeutic opportunity in ovarian cancer* (Bondeson et al.). Each folder is a
self-contained analysis that reads a Figshare dataset, redraws its panels, writes a source-data
workbook, and prints every recomputed number next to the value in the paper (`results/checks_vs_paper.csv`).

| Folder | Panels | Language | Data (Figshare) |
|---|---|---|---|
| [`Figure 1 - Base editor screen`](Figure%201%20-%20Base%20editor%20screen) | Fig 1A-B, ED Fig 1A-E: base editor screen of XPR1 and KIDINS220 | R | https://doi.org/10.6084/m9.figshare.34252485 |
| [`Figure 2 - PRISM chemogenetic profiling`](Figure%202%20-%20PRISM%20chemogenetic%20profiling) | Fig 2D-F, ED Fig 3B-F: PRISM profiling of RBD27 and SC919 | R | https://doi.org/10.6084/m9.figshare.34252608 |
| [`Figure 3 - CellCounter clonality`](Figure%203%20-%20CellCounter%20clonality) | Fig 3D-F, ED Fig 4G: CellCounter clonal barcoding of OVISE xenografts | Python | https://doi.org/10.6084/m9.figshare.34253052 |
| [`Figure 4 - Genome-scale modifier screen`](Figure%204%20-%20Genome-scale%20modifier%20screen) | Fig 4A, ED Fig 5B: genome-scale modifier screens of RBD27 and SC919 | R | https://doi.org/10.6084/m9.figshare.34253274 |

## How to run a figure

1. Download the Figshare dataset of that figure and put its files in the folder's `data/` directory
   (or keep them elsewhere and set `FIG1_DATA_DIR`, `FIG2_DATA_DIR`, `FIG3_DATA_DIR` or `FIG4_DATA_DIR`).
2. Run it as described in the folder's `README.md`: the R analyses are knitted with
   `Rscript -e 'knitr::knit("<file>.Rmd")'`, the Python analysis is `python run_figure3.py`.

The data are not stored in this repository. The `results/` folders hold the expected output.
