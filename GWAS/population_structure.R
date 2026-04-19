#!/usr/bin/env Rscript
#
# Population Structure Analysis for Kelp GWAS Pipeline
#
# This script performs:
# 1. Principal Component Analysis (PCA) with LD pruning
# 2. Kinship/relatedness matrix calculation
# 3. Population differentiation tests (AMOVA, Fst)
# 4. Visualization of population structure
# 5. Check for stratification with phenotype
#
# Author: GWAS Pipeline
# Date: 2026-01-26
#
# Usage:
#   Rscript population_structure.R --plink plink_format \\
#                                  --phenotypes phenotypes.csv \\
#                                  --output-dir population_structure/
#

# Load required libraries
suppressPackageStartupMessages({
  library(SNPRelate)
  library(gdsfmt)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list()
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--plink") {
      params$plink_prefix <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--phenotypes") {
      params$phenotypes <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output-dir") {
      params$output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--n-pcs") {
      params$n_pcs <- as.integer(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}

# Default parameters
params <- list(
  plink_prefix = "plink_format",
  phenotypes = "phenotypes.csv",
  output_dir = "population_structure",
  n_pcs = 10
)

if (length(args) > 0) {
  params <- parse_args(args)
}

# Create output directory
dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)

# Setup logging
log_file <- file.path(params$output_dir, "population_structure.log")
log_conn <- file(log_file, "w")

log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message <- paste0("[", timestamp, "] ", msg)
  cat(message, "\n")
  cat(message, "\n", file = log_conn)
  flush(log_conn)
}

log_message("========================================")
log_message("Population Structure Analysis Pipeline")
log_message("========================================")
log_message(paste("PLINK prefix:", params$plink_prefix))
log_message(paste("Phenotypes:", params$phenotypes))
log_message(paste("Output directory:", params$output_dir))

# ============================================================================
# Step 1: Convert PLINK to GDS format
# ============================================================================
log_message("Step 1: Converting PLINK to GDS format...")

bed_file <- paste0(params$plink_prefix, ".bed")
bim_file <- paste0(params$plink_prefix, ".bim")
fam_file <- paste0(params$plink_prefix, ".fam")

if (!file.exists(bed_file)) {
  stop(paste("BED file not found:", bed_file))
}

gds_file <- file.path(params$output_dir, "genotypes.gds")

tryCatch({
  snpgdsBED2GDS(
    bed.fn = bed_file,
    bim.fn = bim_file,
    fam.fn = fam_file,
    out.gdsfn = gds_file,
    verbose = FALSE
  )
  log_message("GDS file created successfully")
}, error = function(e) {
  log_message(paste("Error converting to GDS:", e$message))
  stop(e)
})

# Open GDS file
genofile <- snpgdsOpen(gds_file)

# ============================================================================
# Step 2: LD Pruning for PCA
# ============================================================================
log_message("Step 2: Performing LD pruning...")

set.seed(42)  # For reproducibility

tryCatch({
  snpset <- snpgdsLDpruning(
    genofile,
    ld.threshold = 0.2,
    slide.max.bp = 50000,
    slide.max.n = 50,
    verbose = FALSE
  )

  # Get pruned SNP list
  pruned_snpset <- unlist(snpset)
  log_message(paste("LD pruning complete:", length(pruned_snpset), "SNPs retained"))

  # Save pruned SNP list
  write.table(
    data.frame(snp_id = pruned_snpset),
    file.path(params$output_dir, "ld_pruned_snps.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

}, error = function(e) {
  log_message(paste("Error in LD pruning:", e$message))
  stop(e)
})

# ============================================================================
# Step 3: Principal Component Analysis
# ============================================================================
log_message("Step 3: Running PCA...")

tryCatch({
  pca_result <- snpgdsPCA(
    genofile,
    snp.id = pruned_snpset,
    num.thread = 4,
    verbose = FALSE
  )

  log_message(paste("PCA complete. Variance explained by first 10 PCs:"))

  # Calculate variance explained
  variance_explained <- pca_result$varprop * 100
  for (i in 1:min(10, length(variance_explained))) {
    log_message(sprintf("  PC%d: %.2f%%", i, variance_explained[i]))
  }

  # Save PC scores
  pc_data <- data.frame(
    sample_id = pca_result$sample.id,
    pca_result$eigenvect[, 1:params$n_pcs]
  )
  colnames(pc_data) <- c("sample_id", paste0("PC", 1:params$n_pcs))

  # Load phenotypes to add population information
  if (file.exists(params$phenotypes)) {
    pheno <- read.csv(params$phenotypes, stringsAsFactors = FALSE)

    # Match samples
    pc_data <- merge(pc_data, pheno, by.x = "sample_id", by.y = "kelp_id", all.x = TRUE)
  }

  write.csv(
    pc_data,
    file.path(params$output_dir, "pca_results.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  # Save variance explained
  variance_df <- data.frame(
    PC = paste0("PC", 1:length(variance_explained)),
    variance_explained = variance_explained
  )
  write.csv(
    variance_df,
    file.path(params$output_dir, "variance_explained.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  log_message("PCA results saved")

}, error = function(e) {
  log_message(paste("Error in PCA:", e$message))
  stop(e)
})

# ============================================================================
# Step 4: Kinship Matrix
# ============================================================================
log_message("Step 4: Calculating kinship matrix...")

tryCatch({
  # Identity-by-state (IBS) matrix
  ibs <- snpgdsIBS(
    genofile,
    snp.id = pruned_snpset,
    num.thread = 4,
    verbose = FALSE
  )

  # Genomic relationship matrix
  grm <- snpgdsGRM(
    genofile,
    snp.id = pruned_snpset,
    num.thread = 4,
    verbose = FALSE
  )

  # Save kinship matrix
  kinship_matrix <- grm$grm
  rownames(kinship_matrix) <- grm$sample.id
  colnames(kinship_matrix) <- grm$sample.id

  write.table(
    kinship_matrix,
    file.path(params$output_dir, "kinship_matrix.txt"),
    quote = FALSE,
    sep = "\t"
  )

  # Save IBS matrix
  ibs_matrix <- ibs$ibs
  rownames(ibs_matrix) <- ibs$sample.id
  colnames(ibs_matrix) <- ibs$sample.id

  write.table(
    ibs_matrix,
    file.path(params$output_dir, "ibs_matrix.txt"),
    quote = FALSE,
    sep = "\t"
  )

  # Identify highly related individuals (kinship > 0.25)
  diag(kinship_matrix) <- 0  # Remove self-relationships
  related_pairs <- which(kinship_matrix > 0.25, arr.ind = TRUE)

  if (nrow(related_pairs) > 0) {
    log_message(paste("Found", nrow(related_pairs), "highly related pairs (kinship > 0.25)"))

    related_df <- data.frame(
      sample1 = rownames(kinship_matrix)[related_pairs[, 1]],
      sample2 = colnames(kinship_matrix)[related_pairs[, 2]],
      kinship = kinship_matrix[related_pairs]
    )

    write.csv(
      related_df,
      file.path(params$output_dir, "related_individuals.csv"),
      row.names = FALSE,
      quote = FALSE
    )
  } else {
    log_message("No highly related individuals found")
  }

  log_message("Kinship analysis complete")

}, error = function(e) {
  log_message(paste("Error calculating kinship:", e$message))
  stop(e)
})

# ============================================================================
# Step 5: Population Statistics
# ============================================================================
log_message("Step 5: Calculating population statistics...")

if ("population" %in% colnames(pc_data)) {
  tryCatch({
    # Calculate Fst between populations
    populations <- unique(pc_data$population[!is.na(pc_data$population)])

    if (length(populations) > 1) {
      log_message("Calculating pairwise Fst...")

      fst_matrix <- matrix(0, nrow = length(populations), ncol = length(populations))
      rownames(fst_matrix) <- populations
      colnames(fst_matrix) <- populations

      for (i in 1:(length(populations) - 1)) {
        for (j in (i + 1):length(populations)) {
          pop1 <- populations[i]
          pop2 <- populations[j]

          # Get sample IDs for each population
          samples_pop1 <- pc_data$sample_id[pc_data$population == pop1 & !is.na(pc_data$population)]
          samples_pop2 <- pc_data$sample_id[pc_data$population == pop2 & !is.na(pc_data$population)]

          # Calculate Fst
          fst <- snpgdsFst(
            genofile,
            population = rep(c(1, 2), c(length(samples_pop1), length(samples_pop2))),
            sample.id = c(samples_pop1, samples_pop2),
            snp.id = pruned_snpset,
            method = "W&C84",
            verbose = FALSE
          )

          fst_matrix[i, j] <- fst$Fst
          fst_matrix[j, i] <- fst$Fst
        }
      }

      write.csv(
        fst_matrix,
        file.path(params$output_dir, "fst_matrix.csv"),
        quote = FALSE
      )

      log_message("Fst calculations complete")
    }
  }, error = function(e) {
    log_message(paste("Warning: Could not calculate Fst:", e$message))
  })
}

# ============================================================================
# Step 6: Check for Stratification with Phenotype
# ============================================================================
log_message("Step 6: Checking for population stratification...")

if ("survival" %in% colnames(pc_data)) {
  tryCatch({
    # Test correlation between PCs and phenotype
    cor_results <- data.frame(
      PC = character(),
      correlation = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )

    for (i in 1:min(10, params$n_pcs)) {
      pc_name <- paste0("PC", i)
      if (pc_name %in% colnames(pc_data)) {
        # Remove NAs
        valid_data <- pc_data[!is.na(pc_data[[pc_name]]) & !is.na(pc_data$survival), ]

        if (nrow(valid_data) > 0) {
          test_result <- cor.test(
            valid_data[[pc_name]],
            valid_data$survival,
            method = "spearman"
          )

          cor_results <- rbind(cor_results, data.frame(
            PC = pc_name,
            correlation = test_result$estimate,
            p_value = test_result$p.value
          ))
        }
      }
    }

    write.csv(
      cor_results,
      file.path(params$output_dir, "pc_phenotype_correlation.csv"),
      row.names = FALSE,
      quote = FALSE
    )

    # Warning if strong correlation
    strong_cor <- cor_results[abs(cor_results$correlation) > 0.3 & cor_results$p_value < 0.05, ]
    if (nrow(strong_cor) > 0) {
      log_message("WARNING: Strong correlation detected between PCs and phenotype!")
      for (i in 1:nrow(strong_cor)) {
        log_message(sprintf("  %s: r=%.3f, p=%.2e",
                           strong_cor$PC[i],
                           strong_cor$correlation[i],
                           strong_cor$p_value[i]))
      }
      log_message("Recommendation: Include these PCs as covariates in GWAS model")
    } else {
      log_message("No strong stratification detected")
    }

    # Check survival difference between populations
    if ("population" %in% colnames(pc_data)) {
      pop_survival <- aggregate(
        survival ~ population,
        data = pc_data[!is.na(pc_data$population) & !is.na(pc_data$survival), ],
        FUN = function(x) c(mean = mean(x), n = length(x))
      )

      log_message("Survival rate by population:")
      for (i in 1:nrow(pop_survival)) {
        log_message(sprintf("  %s: %.1f%% (n=%d)",
                           pop_survival$population[i],
                           pop_survival$survival[i, 1] * 100,
                           pop_survival$survival[i, 2]))
      }
    }

  }, error = function(e) {
    log_message(paste("Warning: Could not test stratification:", e$message))
  })
}

# ============================================================================
# Step 7: Create Visualizations
# ============================================================================
log_message("Step 7: Creating visualizations...")

pdf(file.path(params$output_dir, "pc_plots.pdf"), width = 12, height = 10)

# Plot 1: PC1 vs PC2
p1 <- ggplot(pc_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = if ("population" %in% colnames(pc_data)) population else "All"),
             size = 3, alpha = 0.7) +
  labs(
    title = "Principal Component Analysis",
    subtitle = sprintf("PC1 explains %.2f%%, PC2 explains %.2f%% of variance",
                      variance_explained[1], variance_explained[2]),
    x = sprintf("PC1 (%.2f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.2f%%)", variance_explained[2]),
    color = "Population"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )
print(p1)

# Plot 2: PC1 vs PC3
p2 <- ggplot(pc_data, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = if ("population" %in% colnames(pc_data)) population else "All"),
             size = 3, alpha = 0.7) +
  labs(
    title = "PC1 vs PC3",
    x = sprintf("PC1 (%.2f%%)", variance_explained[1]),
    y = sprintf("PC3 (%.2f%%)", variance_explained[3]),
    color = "Population"
  ) +
  theme_minimal()
print(p2)

# Plot 3: PC2 vs PC3
p3 <- ggplot(pc_data, aes(x = PC2, y = PC3)) +
  geom_point(aes(color = if ("population" %in% colnames(pc_data)) population else "All"),
             size = 3, alpha = 0.7) +
  labs(
    title = "PC2 vs PC3",
    x = sprintf("PC2 (%.2f%%)", variance_explained[2]),
    y = sprintf("PC3 (%.2f%%)", variance_explained[3]),
    color = "Population"
  ) +
  theme_minimal()
print(p3)

# Plot 4: Scree plot
scree_data <- data.frame(
  PC = 1:min(20, length(variance_explained)),
  Variance = variance_explained[1:min(20, length(variance_explained))]
)

p4 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = scree_data$PC)
print(p4)

dev.off()

log_message("PCA plots saved to pc_plots.pdf")

# Kinship heatmap
pdf(file.path(params$output_dir, "kinship_heatmap.pdf"), width = 10, height = 10)

# Sample subset for clearer visualization if too many samples
if (nrow(kinship_matrix) > 100) {
  log_message("Plotting kinship heatmap for subset of samples (too many for full display)")
  sample_indices <- sort(sample(1:nrow(kinship_matrix), 100))
  kinship_subset <- kinship_matrix[sample_indices, sample_indices]
} else {
  kinship_subset <- kinship_matrix
}

pheatmap(
  kinship_subset,
  color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
  main = "Genomic Kinship Matrix",
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 10
)

dev.off()

log_message("Kinship heatmap saved to kinship_heatmap.pdf")

# ============================================================================
# Step 8: Generate Summary Report
# ============================================================================
log_message("Step 8: Generating summary report...")

report_file <- file.path(params$output_dir, "population_stats.txt")
sink(report_file)

cat("========================================\n")
cat("POPULATION STRUCTURE ANALYSIS SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Total samples: %d\n", nrow(pc_data)))
cat(sprintf("Total SNPs: %d\n", length(pruned_snpset)))
cat(sprintf("Number of PCs calculated: %d\n\n", params$n_pcs))

cat("Variance explained by top 10 PCs:\n")
for (i in 1:min(10, length(variance_explained))) {
  cat(sprintf("  PC%d: %.2f%%\n", i, variance_explained[i]))
}
cat(sprintf("  Cumulative (PC1-10): %.2f%%\n\n",
           sum(variance_explained[1:min(10, length(variance_explained))])))

cat("Recommendation for GWAS:\n")
if (variance_explained[1] > 30) {
  cat("  WARNING: PC1 explains >30% of variance - strong population structure\n")
  cat("  RECOMMEND: Include at least 3-5 PCs as covariates in GWAS\n")
} else {
  cat("  Include PC1-PC3 as covariates in GWAS model\n")
}

sink()

log_message(paste("Summary report saved to", report_file))

# Close GDS file
closefn.gds(genofile)

# Close log
log_message("========================================")
log_message("Population Structure Analysis Complete!")
log_message("========================================")
close(log_conn)

cat("\nPopulation structure analysis completed successfully!\n")
cat(paste("Results saved to:", params$output_dir, "\n"))
