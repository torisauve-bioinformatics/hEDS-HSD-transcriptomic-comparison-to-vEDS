**hEDS/HSD RNA-seq Pipeline**

A containerized RNA-seq workflow for processing AmpliSeq-based transcriptomic data, comparing hypermobile Ehlers-Danlos Syndrome (hEDS), Hypermobility Spectrum Disorder (HSD), and control samples. Built with Nextflow and Singularity for reproducible execution on SLURM-managed HPC clusters.

**Overview**

The pipeline has two stages:

1. Nextflow workflow (main.nf): raw FASTQ → trimmed reads → aligned BAM → gene-level counts, with every step running in a pinned Singularity container.
2. R analysis: count matrices → normalized expression → differential expression (limma-voom) → GO enrichment, with PCA and clustering for sample QC.

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

**Data source**

Raw FASTQ files and sample metadata were retrieved using nf-core/fetchngs, which pulls public sequencing data (e.g., from SRA/ENA/GEO) along with a standardized sample sheet. The R analysis expects this sample sheet format (run_accession, fastq_2, scientific_name columns) when parsing sample metadata.

**Requirements**

Nextflow (DSL2)
Singularity (containers pulled automatically from biocontainers — no manual installs needed)
Access to a SLURM cluster (or adjust nextflow.config for a different executor)
R (≥4.x) with the following packages for the downstream analysis: tidyverse, limma, edgeR, pheatmap, clusterProfiler, org.Hs.eg.db

**Pipeline stages**

StepToolContainerWhat it doesTrimcutadaptquay.io/biocontainers/cutadapt:5.2Adapter trimming, minimum length filterAlignSTARquay.io/biocontainers/star:2.7.11bLocal alignment to hg19, sorted BAM outputCountfeatureCounts (Subread)quay.io/biocontainers/subread:2.1.1Gene-level read counting

Resource requests (CPU/memory/time) and SLURM settings are defined per-process in nextflow.config.


**Usage**

bashnextflow run main.nf \
  --ref_dir /path/to/reference \
  --data_dir /path/to/run_directory \
  --slurm_account YOUR_SLURM_ACCOUNT \
  -profile singularity
