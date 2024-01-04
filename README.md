## Introduction

DNA methylation is a crucial epigenetic modification that plays a significant role in regulating gene expression. This repository provides R scripts for the analysis of DNA methylation data generated from the **Infinium HumanMethylationEPIC** and **Infinium Mouse Methylation** arrays.

## Content

### Human EPIC Array Methylation Data

The `human_epic_analysis.R` script is designed to analyze DNA methylation data from the HumanMethylationEPIC BeadChip array.

### Mouse Methylation Array Data

The `mouse_analysis.R` script is tailored for the analysis of DNA methylation data from mouse methylation array experiments.

## Annotation Files

This repository includes annotation files created to analyze the mouse methylation data. These files consist of the manifest file and two additional files providing essential information about genomic features. Find all these files in the `annotation_files` directory.

> **Note:** The manifest file plays a crucial role in the analysis, especially since `minfi` is primarily designed for human data.
