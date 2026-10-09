# Input data

The files in this folder are not stored in the repository. Download the Figshare
dataset **2026 IP6K paper - Figure 1 base editor screen**
(<https://doi.org/10.6084/m9.figshare.34252485>) and place its contents here:

```
base_editor_library_sgRNA_design_GPP.csv
base_editor_library_gene_annotations.csv
poolq_outputs/
    counts-NP_GPP3154_2902200_Bondeson_20211217_p1.txt        raw counts, C>T editor (BE3.9)
    counts-NP_GPP3155_2902200_Bondeson_20211217_p2.txt        raw counts, A>G editor (ABE8e)
    lognorm-NP_GPP3154_..._p1.txt, lognorm-NP_GPP3155_..._p2.txt
    barcodecounts-..., quality-..., correlation-..., runinfo-...   (both plates)
```

`checksums.sha256` in the Figshare dataset lists the SHA-256 of every file. To keep the
files elsewhere, set the environment variable `FIG1_DATA_DIR` to that folder.
