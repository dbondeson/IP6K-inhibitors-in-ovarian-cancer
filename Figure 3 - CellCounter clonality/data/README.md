# Input data

The files in this folder are not stored in the repository. Download the Figshare dataset
**2026 IP6K Paper - Figure 3 CellCounter clonal barcoding** (<https://doi.org/10.6084/m9.figshare.34253052>)
and place its files here:

```
OVISE_xenograft_CellCounter_UMI_read_counts.csv.gz                CellCounter UMI read counts (48 MB)
OVISE_IP_xenograft_animals_randomization_and_outcome.csv          the 15 randomized IP mice and their outcome
OVISE_IP_xenograft_ascites_ex_vivo_BLI.csv                        ex vivo ascites BLI at the end of the study
```

The Figshare README holds a description of every column, and `checksums.sha256` the
SHA-256 of each file. To keep the files elsewhere, set the environment variable
`FIG3_DATA_DIR` to that folder or pass `--data-dir`.
