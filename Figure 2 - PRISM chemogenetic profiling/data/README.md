# Input data

The files in this folder are not stored in the repository. Download the Figshare dataset
**2026 IP6K paper - Figure 2 PRISM screen** (<https://doi.org/10.6084/m9.figshare.34252608>) and place its contents here:

```
sushi_inputs/
    sample_meta.csv, CB_meta.csv, cell_set_meta.csv, cell_line_meta.csv, config.json
sushi_outputs/
    raw_counts.csv.gz               raw barcode counts of every PCR well (not read by the analysis)
    filtered_counts.csv             counts of the PRISM barcodes with sample metadata (ED Fig 3C)
    normalized_counts.csv           counts normalized to the control barcodes (check of the fold changes)
    l2fc.csv                        log2 fold change per condition (check of the fold changes)
    collapsed_l2fc.csv              log2 fold change per cell line, drug and dose (the input of the figures)
PRISM_cell_line_annotations.csv
DepMap_25Q3_mRNA_expression_log2TPMp1_PRISM_cell_lines.csv.gz
DepMap_25Q3_CRISPR_Chronos_gene_effect_PRISM_cell_lines.csv.gz
DepMap_25Q3_CRISPR_Chronos_pathway_genes_all_models.csv
PRISM_AUC_biomarker_correlations_DepMap25Q3_and_repurposing.csv
```

`checksums.sha256` in the Figshare dataset lists the SHA-256 of every file. To keep the files
elsewhere, set the environment variable `FIG2_DATA_DIR` to that folder.
