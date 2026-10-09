# CellCounter clonality of SC919-treated OVISE xenografts

Code for Figure 3D-F and Extended Data Figure 4G of *Dynamic inositol pyrophosphate
synthesis is a targetable therapeutic opportunity in ovarian cancer* (Bondeson et al.).

`run_figure3.py` reads the Figshare dataset for these panels, regenerates them, writes a
source-data workbook, and compares each number quoted in the paper with the value recomputed
from the data. The input is a table of UMI read counts: one row per CellCounter clone (barcode)
per sample, for the 12 mice analysed at the end of the study (5 organs each) and the 15
pre-injection samples.

> **Note on "unique clones" in this submission.** The x axis of Fig 3F and the clone counts in
> the Results (about 26,000 in vehicle and 2,300 in SC919 animals, 11-fold) are the number of
> UMIs detected in each sample (the `n_umis` column), counted before the 3-read filter. The
> CellCounter read threshold is 3 reads (the Methods' "fewer than 10" is an error and will be
> corrected). **Upon re-submission we will recount unique clones from the clones with 3 or more
> reads**, which gives medians of 13,748 (vehicle) and 1,038 (SC919), a 13.2-fold difference.
> Fig 3F, the Results text and the estimate of about 600,000 UMIs in the pre-injection
> samples (about 300,000 with 3 or more reads) will change accordingly.
>
> The current panels and numbers are what the submitted manuscript shows. The workbook sheet
> `unique_clone_sensitivity` gives the medians under each definition.

| Output in `results/` | Panel |
|---|---|
| `fig3d_ascites_BLI_vs_CellCounter_reads.pdf` | Fig 3D: ascites BLI against CellCounter reads |
| `fig3e_clone_size_by_rank.pdf` | Fig 3E: clone size by rank |
| `fig3f_unique_clones_vs_evenness.pdf` | Fig 3F: unique clones and Shannon evenness |
| `ed4g_overlap_treated_lower.pdf`, `ed4g_overlap_preinjection_upper.pdf` | ED Fig 4G: clone overlap between samples (two halves of the panel) |
| `Figure3_CellCounter_source_data.xlsx` | Data behind every panel (sheet guide inside) |
| `checks_vs_paper.csv` | Every recomputed number next to its published value |

The other panels of Figure 3 and Extended Data Figure 4 (pharmacokinetics, whole-body and
organ bioluminescence, tolerability) are not generated here.

## Get the data

Download the Figshare dataset *2026 IP6K Paper - Figure 3 CellCounter clonal barcoding*
(<https://doi.org/10.6084/m9.figshare.34253052>) into `data/`; see `data/README.md` for the file list. To keep it elsewhere:

```bash
export FIG3_DATA_DIR=/path/to/figshare/folder
```

## Run it

```bash
pip install -r requirements.txt
python run_figure3.py
```

Runtime is about 40 seconds and peak memory about 2 GB (4.9 million clone rows). Use `--data-dir`
and `--out-dir` to change the folders. The committed files in `results/` are the expected output.

Tested with Python 3.12.7, pandas 2.2.3, numpy 2.2.5, scipy 1.15.3, matplotlib 3.10.3,
seaborn 0.13.2 and openpyxl 3.2.0b1 on macOS (arm64).

## Analysis

- **Count table.** UMIs were counted from short-read sequencing of the amplified barcode region.
  Spike-in control barcodes were removed before processing, and clones with fewer than 3 reads
  are not in the table. Total reads (Fig 3D), the rank plot (3E), Shannon evenness (3F, y axis)
  and the overlap analysis (ED Fig 4G) all use the clones in the table.
- **Fig 3D.** Spearman correlation of the ex vivo ascites BLI with the summed reads of the
  ascites clones, 12 mice.
- **Fig 3E.** Reads of each clone against its rank within the sample, for the ascites of the 12
  mice and the 15 pre-injection samples.
- **Fig 3F.** *Unique clones* (x axis, and the median clone counts in the Results) is the `n_umis`
  column: the number of UMIs detected in the sample before the 3-read filter (to be replaced by the
  count of clones with 3 or more reads on re-submission, see the note at the top).
  Shannon evenness (y axis) is H' / ln(S) with S the number of clones in the table.
- **ED Fig 4G.** The overlap of two samples is |A ∩ B| / min(|A|, |B|), where A and B are the
  top 10% of the clones (by reads) of each sample. The labels SC1-SC5 and V1-V7 are the mice
  `sH1, sH2, sH5, sI1, sI2` (SC919) and `sJ1-sJ5, sK1, sK2` (vehicle).
- **Mice.** 15 mice were randomized on day 14 (8 SC919, 7 vehicle). Three SC919 mice died in the
  first days of dosing, which leaves 5 SC919 and 7 vehicle mice analysed at the end of the study.

## Checks

`results/checks_vs_paper.csv`: the 9 numbers it lists (randomization and analysis counts, the
Spearman correlation and p-value of Fig 3D, the median unique clones and their fold difference in the
Results, and the sample counts of ED Fig 4G) all reproduce. The p-value of Fig 3D is 8.24e-3, as on the panel;
the legend prints 6.24e-3.

## Layout

```
run_figure3.py     the whole analysis
requirements.txt
data/              put the Figshare files here (README.md lists them)
results/           figures, source-data workbook, checks (expected output)
```
