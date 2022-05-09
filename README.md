## Description

## Content

- [R/MethylationEPIC.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MethylationEPIC.R): Workflow to analyze human methylation array data from the Illumina Infinium Methylation EPIC BeadChip.
- [R/MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R): Script to generate the manifest and annotation files required in [MouseMethylation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseMethylation.R).
- [R/MouseMethylation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseMethylation.R): Workflow to analyze mouse methylation array data from the Illumina Infinium Mouse Methylation BeadChip.
<hasta aquí> 
- [R/plotDMR.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/plotDMR.R): Script to represent the obtained DMRs with a heatmap of the methylation values and their genomic location.

- [data/ann_NCBI_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/ann_NCBI_mouse.csv): File generated with [MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R) including the NCBI annotation of each probe present in the Mouse Methylation array.
- [data/ann_UCSC_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/ann_UCSC_mouse.csv): File generated with [MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R) including the UCSC annotation of each probe present in the Mouse Methylation array.
- [data/enhancers_human.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/enhancers_human.csv) and [data/enhancers_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/enhancers_mouse.csv): Files including FANTOM5 enhancers annotated for human and mouse.
- [data/manifest_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/manifest_mouse.csv): Manifest file for the Mouse Methylation array that maps the probes with their genomic position. It was generated with [MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R).

- [docs/MethylationEPIC.Rmd](https://github.com/raulsanzr/methylation-analysis/blob/main/docs/MethylationEPIC.Rmd): Markdown version of [MethylationEPIC.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MethylationEPIC.R).
- [docs/refs/](https://github.com/raulsanzr/methylation-analysis/blob/main/docs/refs): Folder including the referenced figures.
