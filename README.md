# exorien-exogieo [![DOI](https://zenodo.org/badge/1044976638.svg)](https://doi.org/10.5281/zenodo.16951170)

This repo contains analyses, data, figures and tables relating to the manuscript:

> Jin, N., Hoyd, R., Yilmaz, A. S., Zhu, J., Liu, Y., Singh, M. J., Grencewicz, D., Mo, X., Kalady, M., Rosenberg, D., Dravillas, C. E., Singer, E. A., Carpten, J. D., Chan, C. H., Churchman, M. L., Denko, N. C., Di Clemente, F., Dodd, R. D., Eljilany, I., … Spakowicz, D. (2025). Epigenetic Modulation, Intra-tumoral Microbiome and Immunity in Early Onset Colorectal Cancer. BioRxiv, 2025.03.28.645992. https://doi.org/10.1101/2025.03.28.645992

For example, the script to regenerate the Supplementary Tables is `scripts/tables/tables_supplement.Rmd`, which writes the `.csv` and `.xlsx` files in `tables/`.

```
exorien-exogieo
├── data
├── exorien-exogieo.Rproj
├── figures
├── LICENSE
├── README.md
├── scripts
│   └── tables_supplement.Rmd
└── tables
    ├── compare-deconv-results.csv
    ├── s10_correlation-summary.xlsx
    ├── s6_fig3-rnaseq-vol.xlsx
    ├── s7_fig-16s-vol.xlsx
    ├── s8_orien-immune-frac.csv
    ├── s9_tcga-immune-frac.csv
    ├── tab1_ORIEN_condensed.csv
    ├── tab1_ORIEN.csv
    └── tab1_TCGA_condensed.csv
```

The deconvolved immune cell abundances for the TCGA and ORIEN datasets can be loaded with 

```
loadORIENImmune()
loadTCGAimmune()
```
