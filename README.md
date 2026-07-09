# hEDS/HSD RNA-seq Pipeline

A containerized RNA-seq workflow for processing AmpliSeq-based transcriptomic data, comparing hypermobile Ehlers-Danlos Syndrome (hEDS), Hypermobility Spectrum Disorder (HSD), and control samples. Built with Nextflow and Singularity for reproducible execution on SLURM-managed HPC clusters.

## Overview

The pipeline has two stages:

1. **Nextflow workflow (main.nf):** raw FASTQ → trimmed reads → aligned BAM → gene-level counts, with every step running in a pinned Singularity container.
2. **R analysis:** count matrices → normalized expression → differential expression (limma-voom) → GO enrichment, with PCA and clustering for sample QC.

```
FASTQ files
    │
    ▼
┌─────────┐     cutadapt
│  TRIM   │────────────────► trimmed FASTQ
└─────────┘

    │
    ▼
┌─────────┐     STAR
│  ALIGN  │────────────────► sorted BAM
└─────────┘

    │
    ▼
┌─────────┐     featureCounts
│  COUNT  │────────────────► gene count tables
└─────────┘

    │
    ▼
  R analysis: merge counts → collapse replicates → voom/limma DE → GO enrichment
```


## Data source

Raw FASTQ files and sample metadata were retrieved using nf-core/fetchngs, which pulls public sequencing data (e.g., from SRA/ENA/GEO) along with a standardized sample sheet. The R analysis expects this sample sheet format (run_accession, fastq_2, scientific_name columns) when parsing sample metadata.

## Requirements

Nextflow (DSL2)
Singularity (containers pulled automatically from biocontainers, so no manual installs needed)
Access to a SLURM cluster (or adjust nextflow.config for a different executor)
R (≥4.x) with the following packages for the downstream analysis: tidyverse, limma, edgeR, pheatmap, clusterProfiler, org.Hs.eg.db

## Pipeline stages

| Step | Tool | Container | What it does |
|------|------|-----------|---------------|
| Trim | cutadapt | `quay.io/biocontainers/cutadapt:5.2` | Adapter trimming, minimum length filter |
| Align | STAR | `quay.io/biocontainers/star:2.7.11b` | Local alignment to hg19, sorted BAM output |
| Count | featureCounts (Subread) | `quay.io/biocontainers/subread:2.1.1` | Gene-level read counting |

Resource requests (CPU/memory/time) and SLURM settings are defined per-process in `nextflow.config`.

## Usage

```
bashnextflow run main.nf \
  --ref_dir /path/to/reference \
  --data_dir /path/to/run_directory \
  --slurm_account YOUR_SLURM_ACCOUNT \
  -profile singularity
```
**Required parameters:**
- `--ref_dir` — directory containing `hg19_star_index/` and `hg19.refGene.gtf`
- `--data_dir` — directory containing input FASTQ files; results are written to `data_dir/results`
- `--slurm_account` — SLURM allocation to bill jobs against

Outputs (trimmed reads, BAMs, count tables) are published under `<data_dir>/results/`.

## Downstream analysis (R)

After the Nextflow pipeline produces per-sample count files:

1. **Merge & clean**: count files are combined into a single matrix and matched to sample metadata (condition, patient ID, replicate type) parsed from the nf-core/fetchngs sample sheet.
2. **Replicate handling:** technical replicates are summed; biological replicates are kept separate and modeled with `duplicateCorrelation()`.
3. **Normalization & DE:**: TMM normalization (edgeR) followed by voom/limma with a blocked design, testing hEDS vs. control, HSD vs. control, and hEDS vs. HSD.
4. **QC visualization:** PCA on the top 500 most variable genes, plus k-means clustering to check whether unsupervised groupings align with diagnosis labels.
5. **Functional enrichment:** significant DEGs (FDR < 0.05, |log2FC| > 1) are run through `clusterProfiler::enrichGO()`, with results exported as CSVs and summarized in bar/pie charts by GO category.

Outputs from this stage: `genecounts.csv`, per-comparison `DEGs_*.csv` and `GO_terms_*.csv` files, plus PCA, volcano, and GO summary plots.

## Notes

- Reference genome build: hg19 (RefGene annotation)
- This repo demonstrates the pipeline structure and analysis approach; sample sheets and raw data are not included
- Adjust `nextflow.config` if running on a non-SLURM scheduler or without cluster access (e.g., swap the `executor` block for `local`)

