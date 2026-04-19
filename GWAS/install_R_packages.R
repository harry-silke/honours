#!/usr/bin/env Rscript
#
# R Package Installation Script for Kelp GWAS Pipeline
#
# This script installs all required R packages for the GWAS analysis
#
# Usage:
#   Rscript install_R_packages.R
#
# Author: GWAS Pipeline
# Date: 2026-01-26

cat("=============================================================================\n")
cat("Installing R Packages for Kelp GWAS Pipeline\n")
cat("=============================================================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to install package if not already installed
install_if_missing <- function(pkg, source = "CRAN") {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from", source, "...\n"))

    if (source == "CRAN") {
      install.packages(pkg, dependencies = TRUE)
    } else if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE)
    }

    # Check installation
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("✓", pkg, "installed successfully\n\n"))
    } else {
      cat(paste("✗ Failed to install", pkg, "\n\n"))
    }
  } else {
    cat(paste("✓", pkg, "already installed\n"))
  }
}

# =============================================================================
# Core R Packages (CRAN)
# =============================================================================

cat("\n--- Installing Core Packages ---\n")

cran_packages <- c(
  "data.table",      # Fast data manipulation
  "dplyr",           # Data wrangling
  "tidyr",           # Data tidying
  "ggplot2",         # Visualization
  "reshape2",        # Data reshaping
  "RColorBrewer",    # Color palettes
  "pheatmap",        # Heatmaps
  "gridExtra"        # Arrange multiple plots
)

for (pkg in cran_packages) {
  install_if_missing(pkg, "CRAN")
}

# =============================================================================
# Statistical and Genomics Packages
# =============================================================================

cat("\n--- Installing Statistical Packages ---\n")

stat_packages <- c(
  "MASS",            # Statistical functions
  "lme4",            # Mixed models
  "glmnet",          # Regularized regression
  "pROC",            # ROC curve analysis
  "caret"            # Machine learning
)

for (pkg in stat_packages) {
  install_if_missing(pkg, "CRAN")
}

# =============================================================================
# Survival Analysis Packages
# =============================================================================

cat("\n--- Installing Survival Analysis Packages ---\n")

survival_packages <- c(
  "survival",        # Survival analysis
  "survminer"        # Survival visualization
)

for (pkg in survival_packages) {
  install_if_missing(pkg, "CRAN")
}

# =============================================================================
# Genomic Prediction Packages
# =============================================================================

cat("\n--- Installing Genomic Prediction Packages ---\n")

# rrBLUP (Ridge Regression BLUP)
install_if_missing("rrBLUP", "CRAN")

# =============================================================================
# GWAS Visualization Packages
# =============================================================================

cat("\n--- Installing GWAS Packages ---\n")

# qqman for Manhattan and QQ plots
install_if_missing("qqman", "CRAN")

# =============================================================================
# Bioconductor Packages
# =============================================================================

cat("\n--- Installing Bioconductor Packages ---\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "SNPRelate",       # SNP analysis and PCA
  "gdsfmt",          # GDS file format
  "GenomicRanges",   # Genomic intervals (optional)
  "rtracklayer"      # Import/export genomic annotations (optional)
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, "Bioconductor")
}

# =============================================================================
# Verify Installation
# =============================================================================

cat("\n=============================================================================\n")
cat("Verifying Installation\n")
cat("=============================================================================\n\n")

all_packages <- c(cran_packages, stat_packages, survival_packages,
                 "rrBLUP", "qqman", bioc_packages)

installed <- sapply(all_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

cat("\nInstallation Summary:\n")
cat(paste("Successfully installed:", sum(installed), "/", length(installed), "packages\n\n"))

if (sum(!installed) > 0) {
  cat("Failed to install the following packages:\n")
  cat(paste(" - ", names(installed)[!installed], "\n"))
  cat("\nPlease install these manually or check for errors above.\n")
} else {
  cat("All packages installed successfully! ✓\n")
}

# =============================================================================
# Session Info
# =============================================================================

cat("\n=============================================================================\n")
cat("R Session Information\n")
cat("=============================================================================\n\n")

print(sessionInfo())

# =============================================================================
# Save Package Versions
# =============================================================================

version_file <- "R_package_versions.txt"
sink(version_file)

cat("R Package Versions - Kelp GWAS Pipeline\n")
cat("========================================\n\n")
cat(paste("Date:", Sys.Date(), "\n"))
cat(paste("R version:", R.version.string, "\n\n"))

for (pkg in all_packages[installed]) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(sprintf("%-20s %s\n", pkg, version))
  }
}

sink()

cat("\nPackage versions saved to:", version_file, "\n")
cat("\n=============================================================================\n")
cat("Installation Complete!\n")
cat("=============================================================================\n")
