#!/usr/bin/env python3
"""
VCF Filtering and Processing Module for Kelp GWAS Pipeline

This module performs advanced VCF filtering using cyvcf2:
1. Filter by MAF, missing data, HWE, and depth
2. Keep only biallelic SNPs
3. Calculate per-SNP and per-sample statistics
4. Convert to PLINK format
5. Annotate SNPs with nearest genes
6. Generate QC visualizations

Author: GWAS Pipeline
Date: 2026-01-26
Usage:
    python filter_vcf.py --vcf pass_variants.vcf.gz \\
                        --reference Macpyr2_genomic.fasta \\
                        --annotation Macpyr2_annotation.gff \\
                        --output-dir filtered/ \\
                        --maf 0.05 --missing 0.10
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import subprocess
import pandas as pd
import numpy as np
from scipy import stats
from cyvcf2 import VCF, Writer
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns


class VCFFilter:
    """
    Advanced VCF filtering and processing for GWAS analysis.
    """

    def __init__(self, vcf_file: str, reference: str, annotation: str, output_dir: str,
                 maf: float = 0.05, missing: float = 0.10, hwe_pval: float = 1e-6,
                 min_depth: int = 4, max_depth: int = 20):
        """
        Initialize VCF filtering pipeline.

        Args:
            vcf_file: Input VCF file
            reference: Reference genome FASTA
            annotation: GFF annotation file
            output_dir: Output directory
            maf: Minimum minor allele frequency
            missing: Maximum missing data rate
            hwe_pval: Minimum Hardy-Weinberg equilibrium p-value
            min_depth: Minimum mean depth per SNP
            max_depth: Maximum mean depth per SNP
        """
        self.vcf_file = Path(vcf_file)
        self.reference = Path(reference)
        self.annotation = Path(annotation)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Filtering thresholds
        self.maf_threshold = maf
        self.missing_threshold = missing
        self.hwe_threshold = hwe_pval
        self.min_depth = min_depth
        self.max_depth = max_depth

        # Setup logging
        self.log_file = self.output_dir / "filter_vcf.log"
        self.setup_logging()

        # Statistics storage
        self.snp_stats = []
        self.sample_stats = []
        self.filtering_counts = {
            'input': 0,
            'biallelic': 0,
            'maf': 0,
            'missing': 0,
            'hwe': 0,
            'depth': 0,
            'final': 0
        }

    def setup_logging(self):
        """Configure logging."""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def calculate_maf(self, genotypes: np.ndarray) -> float:
        """
        Calculate minor allele frequency from genotypes.

        Args:
            genotypes: Array of genotypes (-1=missing, 0=ref/ref, 1=ref/alt, 2=alt/alt)

        Returns:
            Minor allele frequency
        """
        # Remove missing genotypes
        genotypes = genotypes[genotypes >= 0]

        if len(genotypes) == 0:
            return 0.0

        # Count alleles
        n_alleles = len(genotypes) * 2
        n_alt = np.sum(genotypes)  # Each genotype contributes 0, 1, or 2 alt alleles

        # Allele frequency
        alt_freq = n_alt / n_alleles
        ref_freq = 1.0 - alt_freq

        # MAF is the minimum
        return min(alt_freq, ref_freq)

    def calculate_missing_rate(self, genotypes: np.ndarray) -> float:
        """
        Calculate missing data rate.

        Args:
            genotypes: Array of genotypes

        Returns:
            Missing rate (0-1)
        """
        n_total = len(genotypes)
        n_missing = np.sum(genotypes < 0)
        return n_missing / n_total if n_total > 0 else 1.0

    def calculate_hwe_pval(self, genotypes: np.ndarray) -> float:
        """
        Calculate Hardy-Weinberg equilibrium p-value.

        Args:
            genotypes: Array of genotypes

        Returns:
            HWE p-value (chi-square test)
        """
        # Remove missing
        genotypes = genotypes[genotypes >= 0]

        if len(genotypes) < 10:  # Too few samples
            return 1.0

        # Count genotypes
        n_00 = np.sum(genotypes == 0)  # Homozygous reference
        n_01 = np.sum(genotypes == 1)  # Heterozygous
        n_11 = np.sum(genotypes == 2)  # Homozygous alternate

        n_total = n_00 + n_01 + n_11

        if n_total == 0:
            return 1.0

        # Calculate allele frequencies
        p = (2 * n_00 + n_01) / (2 * n_total)  # Reference allele frequency
        q = 1.0 - p  # Alternate allele frequency

        # Expected genotype counts under HWE
        exp_00 = n_total * p * p
        exp_01 = n_total * 2 * p * q
        exp_11 = n_total * q * q

        # Chi-square test
        if exp_00 > 0 and exp_01 > 0 and exp_11 > 0:
            chi_sq = ((n_00 - exp_00)**2 / exp_00 +
                     (n_01 - exp_01)**2 / exp_01 +
                     (n_11 - exp_11)**2 / exp_11)

            # 1 degree of freedom for HWE test
            p_val = 1.0 - stats.chi2.cdf(chi_sq, df=1)
            return p_val
        else:
            return 1.0

    def calculate_mean_depth(self, variant) -> float:
        """
        Calculate mean depth for a variant.

        Args:
            variant: cyvcf2 Variant object

        Returns:
            Mean depth across samples
        """
        try:
            depths = variant.format('DP')
            if depths is not None:
                valid_depths = depths[depths >= 0]
                return np.mean(valid_depths) if len(valid_depths) > 0 else 0.0
        except:
            pass
        return 0.0

    def filter_vcf(self) -> Path:
        """
        Apply all filters to VCF file.

        Returns:
            Path to filtered VCF file
        """
        self.logger.info("Starting VCF filtering...")
        self.logger.info(f"Filters: MAF>{self.maf_threshold}, "
                        f"Missing<{self.missing_threshold}, "
                        f"HWE p>{self.hwe_threshold}, "
                        f"Depth {self.min_depth}-{self.max_depth}x")

        # Open input VCF
        vcf_in = VCF(str(self.vcf_file))
        samples = vcf_in.samples

        # Output VCF
        output_vcf = self.output_dir / "final_snps.vcf.gz"
        vcf_out = Writer(str(output_vcf), vcf_in)

        self.filtering_counts['input'] = 0
        variants_written = 0

        # Process each variant
        for variant in vcf_in:
            self.filtering_counts['input'] += 1

            # Progress indicator
            if self.filtering_counts['input'] % 10000 == 0:
                self.logger.info(f"Processed {self.filtering_counts['input']:,} variants, "
                               f"retained {variants_written:,}")

            # Filter 1: Biallelic SNPs only
            if not variant.is_snp or len(variant.ALT) > 1:
                continue
            self.filtering_counts['biallelic'] += 1

            # Get genotypes (-1=missing, 0,1,2=genotypes)
            genotypes = variant.gt_types

            # Filter 2: MAF
            maf = self.calculate_maf(genotypes)
            if maf < self.maf_threshold:
                continue
            self.filtering_counts['maf'] += 1

            # Filter 3: Missing data
            missing_rate = self.calculate_missing_rate(genotypes)
            if missing_rate > self.missing_threshold:
                continue
            self.filtering_counts['missing'] += 1

            # Filter 4: HWE
            hwe_pval = self.calculate_hwe_pval(genotypes)
            if hwe_pval < self.hwe_threshold:
                continue
            self.filtering_counts['hwe'] += 1

            # Filter 5: Depth
            mean_depth = self.calculate_mean_depth(variant)
            if mean_depth < self.min_depth or mean_depth > self.max_depth:
                continue
            self.filtering_counts['depth'] += 1

            # Passed all filters - write to output
            vcf_out.write_record(variant)
            variants_written += 1

            # Store statistics
            self.snp_stats.append({
                'scaffold': variant.CHROM,
                'position': variant.POS,
                'ref': variant.REF,
                'alt': variant.ALT[0],
                'maf': maf,
                'missing_rate': missing_rate,
                'hwe_pval': hwe_pval,
                'mean_depth': mean_depth
            })

        vcf_out.close()
        vcf_in.close()

        self.filtering_counts['final'] = variants_written

        self.logger.info(f"Filtering complete: {variants_written:,} SNPs retained from "
                        f"{self.filtering_counts['input']:,} input variants")

        return output_vcf

    def calculate_sample_statistics(self):
        """Calculate per-sample QC statistics."""
        self.logger.info("Calculating per-sample statistics...")

        vcf = VCF(str(self.output_dir / "final_snps.vcf.gz"))
        samples = vcf.samples

        # Initialize counters
        sample_data = {sample: {'n_genotyped': 0, 'n_het': 0, 'n_total': 0}
                      for sample in samples}

        # Process variants
        for variant in vcf:
            genotypes = variant.gt_types

            for i, sample in enumerate(samples):
                gt = genotypes[i]
                sample_data[sample]['n_total'] += 1

                if gt >= 0:  # Not missing
                    sample_data[sample]['n_genotyped'] += 1
                    if gt == 1:  # Heterozygous
                        sample_data[sample]['n_het'] += 1

        vcf.close()

        # Calculate rates
        for sample in samples:
            data = sample_data[sample]
            genotype_rate = data['n_genotyped'] / data['n_total'] if data['n_total'] > 0 else 0
            het_rate = data['n_het'] / data['n_genotyped'] if data['n_genotyped'] > 0 else 0

            self.sample_stats.append({
                'sample': sample,
                'n_variants': data['n_total'],
                'n_genotyped': data['n_genotyped'],
                'genotype_rate': genotype_rate,
                'heterozygosity': het_rate
            })

        # Identify outliers (possible contamination)
        df = pd.DataFrame(self.sample_stats)
        mean_het = df['heterozygosity'].mean()
        std_het = df['heterozygosity'].std()

        df['het_outlier'] = np.abs(df['heterozygosity'] - mean_het) > 3 * std_het

        outliers = df[df['het_outlier']]
        if len(outliers) > 0:
            self.logger.warning(f"Found {len(outliers)} samples with extreme heterozygosity:")
            for _, row in outliers.iterrows():
                self.logger.warning(f"  {row['sample']}: het={row['heterozygosity']:.3f}")

        self.sample_stats = df.to_dict('records')

    def convert_to_plink(self):
        """Convert VCF to PLINK format using plink."""
        self.logger.info("Converting to PLINK format...")

        input_vcf = self.output_dir / "final_snps.vcf.gz"
        output_prefix = self.output_dir / "plink_format"

        cmd = [
            'plink',
            '--vcf', str(input_vcf),
            '--make-bed',
            '--out', str(output_prefix),
            '--allow-extra-chr',  # For non-standard chromosome names
            '--double-id'  # Use sample ID for both family and individual ID
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("PLINK conversion completed")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"PLINK conversion failed: {e.stderr}")
            self.logger.warning("Continuing without PLINK format...")

    def annotate_nearest_genes(self):
        """Annotate SNPs with nearest genes using bedtools."""
        self.logger.info("Annotating SNPs with nearest genes...")

        # Create BED file from SNPs
        snp_bed = self.output_dir / "snps.bed"
        with open(snp_bed, 'w') as f:
            for snp in self.snp_stats:
                # BED format: chrom, start (0-based), end, name
                f.write(f"{snp['scaffold']}\t{snp['position']-1}\t{snp['position']}\t"
                       f"{snp['scaffold']}:{snp['position']}\n")

        # Extract genes from GFF
        genes_bed = self.output_dir / "genes.bed"
        self.logger.info("Extracting genes from GFF...")

        try:
            with open(self.annotation, 'r') as f_in, open(genes_bed, 'w') as f_out:
                for line in f_in:
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue

                    feature_type = parts[2]
                    if feature_type == 'gene':
                        chrom = parts[0]
                        start = int(parts[3]) - 1  # GFF is 1-based, BED is 0-based
                        end = int(parts[4])
                        attributes = parts[8]

                        # Extract gene ID
                        gene_id = ""
                        for attr in attributes.split(';'):
                            if 'ID=' in attr or 'gene_id=' in attr:
                                gene_id = attr.split('=')[1].strip()
                                break

                        f_out.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")

            # Use bedtools closest to find nearest gene
            annotated_file = self.output_dir / "snp_annotations.txt"

            cmd = [
                'bedtools', 'closest',
                '-a', str(snp_bed),
                '-b', str(genes_bed),
                '-d'  # Report distance
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Parse results and add to SNP stats
            annotations = {}
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 8:
                    snp_id = parts[3]  # scaffold:position
                    gene_id = parts[7] if parts[7] != '.' else 'intergenic'
                    distance = int(parts[8]) if len(parts) > 8 else 0

                    annotations[snp_id] = {
                        'gene_id': gene_id,
                        'distance': distance
                    }

            # Add to SNP stats
            for snp in self.snp_stats:
                snp_id = f"{snp['scaffold']}:{snp['position']}"
                if snp_id in annotations:
                    snp['nearest_gene'] = annotations[snp_id]['gene_id']
                    snp['gene_distance'] = annotations[snp_id]['distance']
                else:
                    snp['nearest_gene'] = 'unknown'
                    snp['gene_distance'] = -1

            self.logger.info("SNP annotation completed")

        except Exception as e:
            self.logger.warning(f"Gene annotation failed: {e}")
            # Add empty annotations
            for snp in self.snp_stats:
                snp['nearest_gene'] = 'unknown'
                snp['gene_distance'] = -1

    def save_statistics(self):
        """Save all statistics to files."""
        self.logger.info("Saving statistics...")

        # SNP statistics
        snp_df = pd.DataFrame(self.snp_stats)
        snp_df.to_csv(self.output_dir / "snp_statistics.csv", index=False)

        # Sample statistics
        sample_df = pd.DataFrame(self.sample_stats)
        sample_df.to_csv(self.output_dir / "sample_statistics.csv", index=False)

        # Filtering summary
        summary_file = self.output_dir / "filtering_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("VCF FILTERING SUMMARY\n")
            f.write("="*80 + "\n\n")

            f.write("Filtering cascade:\n")
            f.write(f"  Input variants: {self.filtering_counts['input']:,}\n")
            f.write(f"  Biallelic SNPs: {self.filtering_counts['biallelic']:,} "
                   f"({100*self.filtering_counts['biallelic']/self.filtering_counts['input']:.1f}%)\n")
            f.write(f"  MAF > {self.maf_threshold}: {self.filtering_counts['maf']:,} "
                   f"({100*self.filtering_counts['maf']/self.filtering_counts['input']:.1f}%)\n")
            f.write(f"  Missing < {self.missing_threshold}: {self.filtering_counts['missing']:,} "
                   f"({100*self.filtering_counts['missing']/self.filtering_counts['input']:.1f}%)\n")
            f.write(f"  HWE p > {self.hwe_threshold}: {self.filtering_counts['hwe']:,} "
                   f"({100*self.filtering_counts['hwe']/self.filtering_counts['input']:.1f}%)\n")
            f.write(f"  Depth {self.min_depth}-{self.max_depth}x: {self.filtering_counts['depth']:,} "
                   f"({100*self.filtering_counts['depth']/self.filtering_counts['input']:.1f}%)\n")
            f.write(f"  Final SNPs: {self.filtering_counts['final']:,} "
                   f"({100*self.filtering_counts['final']/self.filtering_counts['input']:.1f}%)\n\n")

            f.write("SNP statistics:\n")
            f.write(f"  Mean MAF: {snp_df['maf'].mean():.3f}\n")
            f.write(f"  Mean missing rate: {snp_df['missing_rate'].mean():.3f}\n")
            f.write(f"  Mean depth: {snp_df['mean_depth'].mean():.1f}x\n\n")

            f.write("Sample statistics:\n")
            f.write(f"  Mean genotyping rate: {sample_df['genotype_rate'].mean():.3f}\n")
            f.write(f"  Mean heterozygosity: {sample_df['heterozygosity'].mean():.3f}\n")

        self.logger.info(f"Statistics saved to {self.output_dir}")

    def create_visualizations(self):
        """Create QC visualizations."""
        self.logger.info("Creating visualizations...")

        snp_df = pd.DataFrame(self.snp_stats)
        sample_df = pd.DataFrame(self.sample_stats)

        # Set style
        sns.set_style("whitegrid")
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        # 1. MAF distribution
        axes[0, 0].hist(snp_df['maf'], bins=50, edgecolor='black')
        axes[0, 0].set_xlabel('Minor Allele Frequency')
        axes[0, 0].set_ylabel('Count')
        axes[0, 0].set_title('MAF Distribution')
        axes[0, 0].axvline(self.maf_threshold, color='red', linestyle='--',
                          label=f'Threshold={self.maf_threshold}')
        axes[0, 0].legend()

        # 2. Missing data distribution
        axes[0, 1].hist(snp_df['missing_rate'], bins=50, edgecolor='black')
        axes[0, 1].set_xlabel('Missing Rate')
        axes[0, 1].set_ylabel('Count')
        axes[0, 1].set_title('Missing Data Distribution')
        axes[0, 1].axvline(self.missing_threshold, color='red', linestyle='--')

        # 3. Depth distribution
        axes[0, 2].hist(snp_df['mean_depth'], bins=50, edgecolor='black')
        axes[0, 2].set_xlabel('Mean Depth')
        axes[0, 2].set_ylabel('Count')
        axes[0, 2].set_title('Depth Distribution')
        axes[0, 2].axvline(self.min_depth, color='red', linestyle='--')
        axes[0, 2].axvline(self.max_depth, color='red', linestyle='--')

        # 4. Sample genotyping rate
        axes[1, 0].bar(range(len(sample_df)), sample_df['genotype_rate'].values)
        axes[1, 0].set_xlabel('Sample')
        axes[1, 0].set_ylabel('Genotyping Rate')
        axes[1, 0].set_title('Genotyping Rate per Sample')
        axes[1, 0].axhline(0.9, color='red', linestyle='--', label='90%')
        axes[1, 0].legend()

        # 5. Sample heterozygosity
        axes[1, 1].hist(sample_df['heterozygosity'], bins=30, edgecolor='black')
        axes[1, 1].set_xlabel('Heterozygosity')
        axes[1, 1].set_ylabel('Count')
        axes[1, 1].set_title('Heterozygosity Distribution')

        # 6. Filtering cascade
        filter_stages = ['Input', 'Biallelic', 'MAF', 'Missing', 'HWE', 'Depth', 'Final']
        filter_values = [
            self.filtering_counts['input'],
            self.filtering_counts['biallelic'],
            self.filtering_counts['maf'],
            self.filtering_counts['missing'],
            self.filtering_counts['hwe'],
            self.filtering_counts['depth'],
            self.filtering_counts['final']
        ]
        axes[1, 2].bar(range(len(filter_stages)), filter_values, edgecolor='black')
        axes[1, 2].set_xticks(range(len(filter_stages)))
        axes[1, 2].set_xticklabels(filter_stages, rotation=45, ha='right')
        axes[1, 2].set_ylabel('Number of SNPs')
        axes[1, 2].set_title('Filtering Cascade')

        plt.tight_layout()
        plt.savefig(self.output_dir / "qc_plots.png", dpi=300, bbox_inches='tight')
        plt.close()

        self.logger.info(f"Visualizations saved to {self.output_dir / 'qc_plots.png'}")

    def run_pipeline(self):
        """Run the complete filtering pipeline."""
        self.logger.info("="*80)
        self.logger.info("Starting VCF Filtering Pipeline")
        self.logger.info(f"Input VCF: {self.vcf_file}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info("="*80)

        # Step 1: Filter VCF
        self.filter_vcf()

        # Step 2: Calculate sample statistics
        self.calculate_sample_statistics()

        # Step 3: Convert to PLINK
        self.convert_to_plink()

        # Step 4: Annotate with nearest genes
        if self.annotation.exists():
            self.annotate_nearest_genes()
        else:
            self.logger.warning(f"Annotation file not found: {self.annotation}")

        # Step 5: Save statistics
        self.save_statistics()

        # Step 6: Create visualizations
        self.create_visualizations()

        self.logger.info("="*80)
        self.logger.info("VCF Filtering Pipeline Completed!")
        self.logger.info(f"Final SNPs: {self.filtering_counts['final']:,}")
        self.logger.info(f"Results saved to: {self.output_dir}")
        self.logger.info("="*80)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='VCF Filtering Pipeline for Kelp GWAS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python filter_vcf.py --vcf pass_variants.vcf.gz \\
                        --reference Macpyr2_genomic.fasta \\
                        --annotation Macpyr2_annotation.gff \\
                        --output-dir filtered/ \\
                        --maf 0.05 --missing 0.10
        """
    )

    parser.add_argument('--vcf', required=True,
                       help='Input VCF file (gzipped)')
    parser.add_argument('--reference', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('--annotation', required=True,
                       help='GFF annotation file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('--maf', type=float, default=0.05,
                       help='Minimum minor allele frequency (default: 0.05)')
    parser.add_argument('--missing', type=float, default=0.10,
                       help='Maximum missing data rate (default: 0.10)')
    parser.add_argument('--hwe-pval', type=float, default=1e-6,
                       help='Minimum HWE p-value (default: 1e-6)')
    parser.add_argument('--min-depth', type=int, default=4,
                       help='Minimum mean depth (default: 4)')
    parser.add_argument('--max-depth', type=int, default=20,
                       help='Maximum mean depth (default: 20)')

    args = parser.parse_args()

    # Create and run pipeline
    pipeline = VCFFilter(
        vcf_file=args.vcf,
        reference=args.reference,
        annotation=args.annotation,
        output_dir=args.output_dir,
        maf=args.maf,
        missing=args.missing,
        hwe_pval=args.hwe_pval,
        min_depth=args.min_depth,
        max_depth=args.max_depth
    )

    pipeline.run_pipeline()


if __name__ == '__main__':
    main()
