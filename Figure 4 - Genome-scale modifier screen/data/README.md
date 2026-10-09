# Input data

The files in this folder are not stored in the repository. Download the Figshare dataset
**2026 IP6K paper - Figure 4 genome-scale drug modifier screen**
(<https://doi.org/10.6084/m9.figshare.34253274>) and place its contents here:

```
RMGI_Humagne_SetC_summed_read_counts.csv                  read counts, library set C (input of the analysis)
RMGI_Humagne_SetD_summed_read_counts.csv                  read counts, library set D (input of the analysis)
RMGI_Humagne_sample_key.csv                               how each count column was assembled
Humagne_construct_to_gene_map.csv                         the 40,376 scored constructs and their genes
GPP_portal_hypergeometric_RBD27_vs_DMSO.txt               GPP portal gene scores (verification only)
GPP_portal_hypergeometric_SC919_vs_DMSO.txt
DepMap_25Q3_RMGI_CRISPR_gene_effect.csv                   RMGI Chronos gene effect, DepMap 25Q3
```

`checksums.sha256` in the Figshare dataset lists the SHA-256 of every file. To keep the files
elsewhere, set the environment variable `FIG4_DATA_DIR` to that folder.
