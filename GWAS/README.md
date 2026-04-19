# Kelp Heat Stress GWAS Pipeline

Comprehensive genome-wide association study (GWAS) pipeline for analyzing heat stress survival in giant kelp (*Macrocystis pyrifera*).

![Pipeline Status](https://img.shields.io/badge/status-production--ready-brightgreen)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![R](https://img.shields.io/badge/R-4.0%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Pipeline Architecture](#pipeline-architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Overview

This pipeline performs a complete GWAS analysis from raw sequencing reads to candidate gene identification and genomic prediction. It is specifically designed for analyzing binary survival phenotypes in kelp after heat stress treatment but can be adapted for other organisms and phenotypes.

**Organism:** *Macrocystis pyrifera* (Giant Kelp)
**Reference Genome:** Macpyr2 from JGI Phycocosm
**Genome Size:** 537 MB, 35 scaffolds, 25,919 genes
**Phenotype:** Binary survival (survived/died after 24°C for 4 weeks)
**Expected Sample Size:** 500 individuals

---

## Features

### Comprehensive Analysis Pipeline
- ✅ Quality control and read alignment
- ✅ GATK-based variant calling
- ✅ Advanced VCF filtering (MAF, HWE, depth, missing data)
- ✅ Population structure analysis (PCA, kinship)
- ✅ Logistic regression GWAS with covariate control
- ✅ Publication-quality visualizations
- ✅ Functional gene annotation
- ✅ Genomic prediction and breeding value estimation

### Robust Design
- ✅ Checkpoint system for resuming interrupted analyses
- ✅ Parallel processing for computational efficiency
- ✅ Comprehensive logging and error handling
- ✅ Configurable parameters via YAML
- ✅ Production-ready code with extensive documentation

### Statistical Rigor
- ✅ Genomic inflation control (lambda calculation)
- ✅ Multiple testing correction (Bonferroni, FDR)
- ✅ Population stratification detection
- ✅ Cross-validated genomic prediction (AUC, accuracy)
- ✅ Effect size estimation with confidence intervals

---

## Pipeline Architecture

The pipeline consists of **9 modules** organized in a sequential workflow:

```
┌─────────────────────────────────────────────────────────────┐
│                    RAW FASTQ FILES                          │
└────────────────────┬────────────────────────────────────────┘
                     │
         ┌───────────▼───────────┐
         │   1. QC & ALIGNMENT   │  FastQC, fastp, HISAT2
         └───────────┬───────────┘
                     │
         ┌───────────▼───────────┐
         │  2. VARIANT CALLING   │  GATK HaplotypeCaller
         └───────────┬───────────┘
                     │
         ┌───────────▼───────────┐
         │   3. VCF FILTERING    │  MAF, HWE, depth, missing
         └───────────┬───────────┘
                     │
      ┌──────────────┴──────────────┐
      │                             │
┌─────▼──────┐              ┌──────▼───────┐
│4. POP      │              │  5. GWAS     │  Logistic regression
│ STRUCTURE  │              │  ANALYSIS    │  + mixed model
└─────┬──────┘              └──────┬───────┘
      │                             │
      │     ┌───────────────────────┘
      │     │
┌─────▼─────▼───┐      ┌──────────────┐      ┌──────────────┐
│6. VISUALIZE   │      │7. ANNOTATE   │      │8. GENOMIC    │
│   GWAS        │      │  GENES       │      │ PREDICTION   │
└───────────────┘      └──────────────┘      └──────────────┘
```

**Python Modules:**
- `qc_alignment.py` - Quality control and alignment
- `variant_calling.py` - GATK variant calling wrapper
- `filter_vcf.py` - VCF filtering with cyvcf2
- `annotate_genes.py` - Gene annotation and enrichment
- `run_gwas_pipeline.py` - Main orchestrator

**R Modules:**
- `population_structure.R` - PCA and kinship analysis
- `gwas_analysis.R` - Association testing
- `visualize_gwas.R` - Manhattan plots, QQ plots
- `genomic_prediction.R` - Prediction models and breeding values

---

## Installation

### Option 1: Conda Installation (Recommended)

```bash
# Install dependencies
bash install_dependencies.sh --method conda

# Activate environment
conda activate kelp_gwas

# Verify installation
bash install_dependencies.sh --method verify
```

### Option 2: Manual Installation

```bash
# Install Python dependencies
pip install -r requirements.txt

# Install R packages
Rscript install_R_packages.R

# Install external tools
# See install_dependencies.sh for tool-specific instructions
```

### System Requirements

- **OS:** Linux, macOS, or Windows WSL2
- **RAM:** 32 GB minimum, 64 GB recommended
- **Disk Space:** 1 TB for 500 samples with ~300 GB sequencing data
- **CPU:** 8+ cores recommended

---

## Quick Start

### 1. Prepare Your Data

```bash
# Create project directory
mkdir kelp_gwas_project
cd kelp_gwas_project

# Copy pipeline scripts
cp /path/to/pipeline/*.py .
cp /path/to/pipeline/*.R .
cp /path/to/pipeline/config.yaml .

# Organize input data
mkdir -p raw_data reference data
# Place FASTQ files in raw_data/
# Place reference genome in reference/
# Place phenotype file in data/
```

### 2. Edit Configuration

Edit `config.yaml` to specify your file paths and parameters:

```yaml
input:
  fastq_dir: "raw_data/"
  reference_genome: "reference/Macpyr2_genomic.fasta"
  annotation: "reference/Macpyr2_annotation.gff"
  phenotypes: "data/phenotypes.csv"

output:
  base_dir: "gwas_output/"

resources:
  threads: 8
  parallel_samples: 4
```

### 3. Run Pipeline

```bash
# Run complete pipeline
python run_gwas_pipeline.py --config config.yaml

# Or run individual modules
python qc_alignment.py --fastq-dir raw_data/ --reference reference/Macpyr2_genomic.fasta --output-dir output/qc/
python variant_calling.py --bam-dir output/aligned/ --reference reference/Macpyr2_genomic.fasta --output-dir output/variants/
# ... etc
```

### 4. View Results

```bash
# Check main output directory
ls gwas_output/

# View significant SNPs
head gwas_output/06_gwas_results/gwas_significant.txt

# View figures
open gwas_output/07_figures/manhattan_plot.pdf

# Read summary report
cat gwas_output/final_report.txt
```

---

## Detailed Usage

### Module 1: Quality Control and Alignment

```bash
python qc_alignment.py \
  --fastq-dir raw_data/ \
  --reference reference/Macpyr2_genomic.fasta \
  --output-dir output/qc/ \
  --threads 8 \
  --parallel 4
```

**Outputs:**
- `qc_reports/` - FastQC HTML reports
- `trimmed/` - Quality-filtered FASTQ files
- `aligned/` - Sorted BAM files with indices
- `qc_summary.csv` - QC metrics for all samples

**Key Parameters:**
- `--threads`: Threads per alignment job
- `--parallel`: Number of samples to process simultaneously

### Module 2: Variant Calling

```bash
python variant_calling.py \
  --bam-dir output/aligned/ \
  --reference reference/Macpyr2_genomic.fasta \
  --output-dir output/variants/ \
  --threads 8 \
  --parallel 10
```

**Outputs:**
- `gvcfs/` - Per-sample GVCF files
- `genomicsdb/` - GenomicsDB workspace
- `vcfs/pass_variants.vcf.gz` - Filtered variants
- `variant_stats.txt` - Variant statistics

### Module 3: VCF Filtering

```bash
python filter_vcf.py \
  --vcf output/variants/vcfs/pass_variants.vcf.gz \
  --reference reference/Macpyr2_genomic.fasta \
  --annotation reference/Macpyr2_annotation.gff \
  --output-dir output/filtered/ \
  --maf 0.05 \
  --missing 0.10
```

**Outputs:**
- `final_snps.vcf.gz` - High-quality SNPs (~600,000 expected)
- `plink_format.*` - PLINK binary format files
- `snp_statistics.csv` - Per-SNP quality metrics
- `sample_statistics.csv` - Per-sample quality metrics

**Filtering Cascade:**
1. Biallelic SNPs only
2. MAF > 0.05 (5%)
3. Missing data < 10%
4. HWE p-value > 1×10⁻⁶
5. Depth 4-20×

### Module 4: Population Structure

```bash
Rscript population_structure.R \
  --plink output/filtered/plink_format \
  --phenotypes data/phenotypes.csv \
  --output-dir output/pop_structure/ \
  --n-pcs 10
```

**Outputs:**
- `pca_results.csv` - PC scores for all samples
- `kinship_matrix.txt` - Genomic relationship matrix
- `pc_plots.pdf` - PCA visualizations
- `variance_explained.csv` - Variance per PC

**Quality Checks:**
- Genomic inflation should be 0.95-1.05
- PC1 should not explain >30% variance
- Check for stratification with phenotype

### Module 5: GWAS Analysis

```bash
Rscript gwas_analysis.R \
  --geno output/filtered/plink_format \
  --pheno data/phenotypes.csv \
  --pcs output/pop_structure/pca_results.csv \
  --kinship output/pop_structure/kinship_matrix.txt \
  --output-dir output/gwas/ \
  --n-pcs 3
```

**Outputs:**
- `gwas_all_snps.txt` - Results for all SNPs
- `gwas_significant.txt` - Genome-wide significant SNPs
- `gwas_top100.txt` - Top 100 SNPs by p-value
- `lambda.txt` - Genomic inflation factor

**Model:** `survival ~ SNP + PC1 + PC2 + PC3 + initial_size_cm`

**Significance Thresholds:**
- Genome-wide: p < 5×10⁻⁸ (Bonferroni)
- Suggestive: p < 1×10⁻⁵

### Module 6: Visualization

```bash
Rscript visualize_gwas.R \
  --gwas output/gwas/gwas_all_snps.txt \
  --pheno data/phenotypes.csv \
  --output figures/
```

**Outputs:**
- `manhattan_plot.pdf` - Genome-wide Manhattan plot
- `qq_plot.pdf` - QQ plot with lambda
- `effect_sizes.pdf` - Forest plot of top SNPs
- `maf_vs_effect.png` - MAF vs significance scatter

### Module 7: Gene Annotation

```bash
python annotate_genes.py \
  --snps output/gwas/gwas_significant.txt \
  --annotation reference/Macpyr2_annotation.gff \
  --output-dir output/annotation/ \
  --window 50000
```

**Outputs:**
- `snp_annotations.csv` - Detailed SNP annotations
- `priority_genes.csv` - Prioritized candidate genes
- `functional_enrichment.txt` - Enrichment test results
- `annotation_plots.png` - Visualization

**Functional Categories:**
- Heat shock & stress response
- Photosynthesis
- Antioxidant defense
- Protein degradation
- Membrane & lipid metabolism

### Module 8: Genomic Prediction

```bash
Rscript genomic_prediction.R \
  --geno output/filtered/plink_format \
  --pheno data/phenotypes.csv \
  --gwas output/gwas/gwas_all_snps.txt \
  --output output/prediction/
```

**Outputs:**
- `model_comparison.csv` - AUC for each model
- `gebvs_all_kelps.csv` - Breeding values for all individuals
- `top_candidates.csv` - Top 20% for breeding/restoration
- `roc_curves.pdf` - Model performance comparison

**Models Tested:**
1. Ridge Regression BLUP
2. Elastic Net
3. Top SNPs only
4. Weighted by effect size

**Expected Performance:** AUC 0.70-0.80

---

## Input Files

### 1. FASTQ Files

Paired-end Illumina reads (150 bp), whole-genome sequencing

```
raw_data/
├── Sample_001_R1.fastq.gz
├── Sample_001_R2.fastq.gz
├── Sample_002_R1.fastq.gz
├── Sample_002_R2.fastq.gz
└── ...
```

### 2. Reference Genome

FASTA format (Macpyr2_genomic.fasta)

### 3. Gene Annotation

GFF format (Macpyr2_annotation.gff)

### 4. Phenotype File

CSV format with required columns:

```csv
kelp_id,survival,initial_size_cm,tank_id,population,week_died
K001,1,25.3,T1,San Diego,NA
K002,0,23.8,T1,San Diego,3
K003,1,28.1,T2,Santa Barbara,NA
...
```

**Required columns:**
- `kelp_id`: Unique sample identifier (must match FASTQ file names)
- `survival`: Binary (1 = survived, 0 = died)

**Optional columns:**
- `initial_size_cm`: Initial frond length
- `tank_id`: Experimental tank
- `population`: Geographic origin
- `week_died`: Week of death (for survival analysis)

---

## Output Files

The pipeline creates a structured output directory:

```
gwas_output/
├── 01_qc/
│   ├── qc_reports/          # FastQC HTML reports
│   ├── trimmed/             # Filtered FASTQ files
│   ├── aligned/             # BAM files
│   └── qc_summary.csv
├── 02_variant_calling/
│   ├── gvcfs/               # Per-sample GVCFs
│   ├── genomicsdb/          # Combined database
│   └── vcfs/                # Joint-called VCF
├── 03_filtering/
│   ├── final_snps.vcf.gz   # ~600K high-quality SNPs
│   ├── plink_format.*       # PLINK files
│   └── snp_statistics.csv
├── 04_population_structure/
│   ├── pca_results.csv
│   ├── kinship_matrix.txt
│   └── pc_plots.pdf
├── 05_gwas_results/
│   ├── gwas_all_snps.txt   # All SNP results
│   ├── gwas_significant.txt # Significant SNPs
│   └── lambda.txt           # Genomic inflation
├── 06_figures/
│   ├── manhattan_plot.pdf
│   ├── qq_plot.pdf
│   ├── effect_sizes.pdf
│   └── ...
├── 07_annotation/
│   ├── snp_annotations.csv
│   ├── priority_genes.csv
│   └── functional_enrichment.txt
├── 08_genomic_prediction/
│   ├── gebvs_all_kelps.csv
│   ├── top_candidates.csv
│   └── model_comparison.csv
└── final_report.html        # Comprehensive HTML report
```

---

## Configuration

Edit `config.yaml` to customize pipeline behavior. Key sections:

### Quality Control
```yaml
qc:
  min_quality: 20
  min_length: 50
  min_mapping_rate: 80
  min_depth: 4
```

### VCF Filtering
```yaml
vcf_filtering:
  min_maf: 0.05
  max_missing: 0.10
  hwe_pval: 0.000001
```

### GWAS Analysis
```yaml
gwas:
  n_pcs: 3
  covariates:
    - initial_size_cm
  use_mixed_model: true
```

### Computational Resources
```yaml
resources:
  threads: 8
  parallel_samples: 4
  memory_gb: 64
```

---

## Troubleshooting

### Common Issues

**1. Out of Memory**
- Reduce `parallel_samples` in config.yaml
- Process samples in batches
- Use a high-memory node if on HPC

**2. GATK Errors**
- Ensure Java 8+ is installed
- Check GATK version (4.0+)
- Verify reference genome has .dict and .fai files

**3. Low Mapping Rate (<80%)**
- Check FASTQ quality with FastQC
- Verify correct reference genome
- Consider more aggressive quality filtering

**4. High Genomic Inflation (λ > 1.1)**
- Include more PCs as covariates (increase `n_pcs`)
- Use mixed model GWAS (set `use_mixed_model: true`)
- Check for cryptic relatedness

**5. No Significant SNPs**
- Check sample size (need ~200+ for 80% power)
- Verify phenotype variance
- Consider less stringent significance threshold
- Check for batch effects

### Getting Help

1. Check log files in each output directory
2. Verify input file formats
3. Test with example data
4. Open an issue on GitHub
5. Contact: hisilke@utas.edu.au

---

## Performance Benchmarks

**Test System:** 16-core CPU, 64 GB RAM, 1 TB SSD

| Step | Time (500 samples) | Disk Space |
|------|-------------------|------------|
| QC & Alignment | 24-48 hours | 200 GB |
| Variant Calling | 48-72 hours | 100 GB |
| Filtering | 2-4 hours | 10 GB |
| GWAS | 1-2 hours | 5 GB |
| Visualization | <1 hour | 1 GB |

**Total Runtime:** ~5-7 days
**Total Disk Space:** ~300-400 GB

---

## Citation

If you use this pipeline in your research, please cite:

```
Harry Silke (2026). Kelp Heat Stress GWAS Pipeline.
GitHub: https://github.com/harry-silke/utas-honours
```

And cite the tools used:
- **GATK:** McKenna et al. (2010) Genome Research
- **PLINK:** Chang et al. (2015) GigaScience
- **SNPRelate:** Zheng et al. (2012) Bioinformatics
- **cyvcf2:** Pedersen & Quinlan (2017) Bioinformatics

---

## License

MIT License - See LICENSE file for details

---

## Authors

**Principal Investigator:** Harry Silke
**Email:** hisilke@utas.edu.au
**Institution:** IMAS
**Date:** January 2026

---

## Acknowledgments

- JGI Phycocosm for *Macrocystis pyrifera* genome
- GATK team at Broad Institute
- SNPRelate and Bioconductor community
- All contributors to open-source bioinformatics tools

---

**Version:** 1.0.0
**Last Updated:** January 26, 2026
**Status:** Production-Ready
