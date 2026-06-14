# Small RNA sequencing analysis of serum and sperm samples associated with adverse childhood experiences

This repository contains R code used for secondary analysis of human small RNA sequencing data from serum and sperm samples.

## Data availability

Raw sequencing files and processed count data are deposited in ArrayExpress/BioStudies under accession number: E-MTAB-17125 https://doi.org/10.6019/E-MTAB-17125.

Sensitive participant-level metadata are not included in this repository.

## Repository contents

- `scripts/01_miRNAseq_secondary_analysis_cleaned.R` — R script for quality assessment, filtering, differential expression analysis, and visualization of small RNA sequencing data.

## Analysis overview

The analysis includes:

- loading primary QIAGEN GeneGlobe UMI count tables,
- removal of low-quality samples,
- filtering of lowly expressed miRNAs,
- differential expression analysis using DESeq2,
- PCA and sample correlation analysis,
- visualization of selected results,
- optional exploratory miRNA target enrichment analysis.

## Required input files

The script expects input files in a local `input/` directory. These files are not included in this repository because they contain sequencing-derived data and/or sensitive participant-level metadata.

Expected local input files:

- `primary_UMI_counts.csv`
- `illumina_sample_registry.tsv`
- `sample_metadata.tsv`
- optionally: `tarbase_data.tsv`

## Output

The script creates:

- `results/` — differential expression tables and session information,
- `figures/` — PCA plots, correlation plots, heatmaps, volcano plots, and normalized count plots.

## Software

The analysis was performed in R using DESeq2, ggplot2, pheatmap, corrplot, EnhancedVolcano, and related Bioconductor packages.

## Notes

This code was adapted from the original analysis script for documentation and reproducibility purposes. File paths and input file names should be adjusted locally before running the script.
