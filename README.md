## Description

## Content

### [/R](https://github.com/raulsanzr/methylation-analysis/tree/main/R)

- [MethylationEPIC.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MethylationEPIC.R): Workflow to analyze human methylation array data from the Infinium MethylationEPIC BeadChip.
- [MouseMethylation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseMethylation.R): Workflow to analyze mouse methylation array data from the Infinium Mouse Methylation BeadChip.
- [MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R): Script to generate the manifest and annotation files required in [MouseMethylation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseMethylation.R).
- [plotDMR.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/plotDMR.R): Script to represent the DMRs obtained with the two first scripts together with a heatmap of the methylation values and its genomic location.

### [/data](https://github.com/raulsanzr/methylation-analysis/tree/main/data)

- [ann_NCBI_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/ann_NCBI_mouse.csv): Annotation file produced with [MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R) including the NCBI annotation of each probe present in the Mouse Methylation array.
- [ann_UCSC_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/ann_UCSC_mouse.csv): Annotation file produced with [MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R) including the UCSC annotation of each probe present in the Mouse Methylation array.
- [manifest_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/manifest_mouse.csv): Manifest file for the Mouse Methylation array mapping the probes with its genomic position produced with [MouseAnnotation.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MouseAnnotation.R).
- [enhancers_human.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/enhancers_human.csv) and [enhancers_mouse.csv](https://github.com/raulsanzr/methylation-analysis/blob/main/data/enhancers_mouse.csv): Files including FANTOM5 enhancers annotated for human and mouse.

### [/docs](https://github.com/raulsanzr/methylation-analysis/tree/main/docs)

- [MethylationEPIC.Rmd](https://github.com/raulsanzr/methylation-analysis/blob/main/docs/MethylationEPIC.Rmd): Markdown version of [MethylationEPIC.R](https://github.com/raulsanzr/methylation-analysis/blob/main/R/MethylationEPIC.R).
- [/refs](https://github.com/raulsanzr/methylation-analysis/blob/main/docs/refs): Folder including the referenced figures in the Rmd files.
