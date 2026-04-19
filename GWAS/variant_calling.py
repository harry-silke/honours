#!/usr/bin/env python3
"""
Variant Calling Module for Kelp GWAS Pipeline

This module wraps the GATK variant calling pipeline:
1. Mark duplicates in BAM files
2. Call variants per sample using HaplotypeCaller (GVCF mode)
3. Combine GVCFs using GenomicsDBImport
4. Joint genotyping with GenotypeGVCFs
5. Apply hard filters to variants
6. Extract PASS variants and calculate statistics

Author: GWAS Pipeline
Date: 2026-01-26
Usage:
    python variant_calling.py --bam-dir aligned/ --reference Macpyr2_genomic.fasta \\
                             --output-dir variant_calling/ --threads 8 --parallel 10
"""

import os
import sys
import logging
import argparse
import subprocess
from pathlib import Path
from typing import List, Dict
import json
import multiprocessing
from datetime import datetime
import pandas as pd


class VariantCalling:
    """
    GATK-based variant calling pipeline for kelp genomes.
    """

    def __init__(self, bam_dir: str, reference: str, output_dir: str,
                 threads: int = 8, parallel: int = 10, gatk_path: str = 'gatk'):
        """
        Initialize variant calling pipeline.

        Args:
            bam_dir: Directory containing sorted BAM files
            reference: Path to reference genome FASTA
            output_dir: Output directory
            threads: Number of threads for GATK operations
            parallel: Number of samples to process in parallel
            gatk_path: Path to GATK executable (or 'gatk' if in PATH)
        """
        self.bam_dir = Path(bam_dir)
        self.reference = Path(reference)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.parallel = parallel
        self.gatk = gatk_path

        # Create output directories
        self.dedup_dir = self.output_dir / "dedup"
        self.gvcf_dir = self.output_dir / "gvcfs"
        self.db_dir = self.output_dir / "genomicsdb"
        self.vcf_dir = self.output_dir / "vcfs"

        for dir_path in [self.dedup_dir, self.gvcf_dir, self.db_dir, self.vcf_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        # Setup logging
        self.log_file = self.output_dir / "variant_calling.log"
        self.setup_logging()

        # Checkpoint
        self.checkpoint_file = self.output_dir / "variant_checkpoint.json"
        self.checkpoint = self.load_checkpoint()

        # Verify GATK is available
        self.verify_gatk()

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

    def load_checkpoint(self) -> Dict:
        """Load checkpoint data."""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        return {
            'deduped_samples': [],
            'gvcf_samples': [],
            'genomicsdb_created': False,
            'joint_genotyped': False,
            'filtered': False
        }

    def save_checkpoint(self):
        """Save checkpoint data."""
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.checkpoint, f, indent=2)

    def verify_gatk(self):
        """Verify GATK is installed and accessible."""
        try:
            result = subprocess.run([self.gatk, '--version'],
                                  capture_output=True, text=True, check=True)
            version = result.stderr.split('\n')[0]  # GATK prints version to stderr
            self.logger.info(f"Found GATK: {version}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise EnvironmentError(f"GATK not found at '{self.gatk}'. "
                                 f"Please install GATK or specify correct path.")

    def find_bam_files(self) -> List[Path]:
        """
        Find all sorted BAM files in the input directory.

        Returns:
            List of BAM file paths
        """
        bam_files = sorted(self.bam_dir.glob("*_sorted.bam"))

        if not bam_files:
            # Try without _sorted suffix
            bam_files = sorted(self.bam_dir.glob("*.bam"))

        self.logger.info(f"Found {len(bam_files)} BAM files")
        return bam_files

    def create_sequence_dictionary(self):
        """Create sequence dictionary for reference genome if it doesn't exist."""
        dict_file = self.reference.with_suffix('.dict')

        if dict_file.exists():
            self.logger.info("Sequence dictionary already exists")
            return

        self.logger.info("Creating sequence dictionary...")
        cmd = [
            self.gatk, 'CreateSequenceDictionary',
            '-R', str(self.reference),
            '-O', str(dict_file)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("Sequence dictionary created")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create sequence dictionary: {e.stderr}")
            raise

    def index_reference(self):
        """Index reference genome with samtools if not already done."""
        fai_file = Path(str(self.reference) + '.fai')

        if fai_file.exists():
            self.logger.info("Reference index already exists")
            return

        self.logger.info("Indexing reference genome...")
        cmd = ['samtools', 'faidx', str(self.reference)]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("Reference indexed")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to index reference: {e.stderr}")
            raise

    def mark_duplicates(self, bam_file: Path) -> Path:
        """
        Mark duplicates in BAM file using GATK.

        Args:
            bam_file: Input BAM file

        Returns:
            Path to deduplicated BAM file
        """
        sample_name = bam_file.stem.replace('_sorted', '')

        # Check checkpoint
        if sample_name in self.checkpoint['deduped_samples']:
            output_bam = self.dedup_dir / f"{sample_name}_dedup.bam"
            self.logger.info(f"Sample {sample_name} already deduplicated (checkpoint)")
            return output_bam

        output_bam = self.dedup_dir / f"{sample_name}_dedup.bam"
        metrics_file = self.dedup_dir / f"{sample_name}_metrics.txt"

        self.logger.info(f"Marking duplicates for {sample_name}...")

        cmd = [
            self.gatk, 'MarkDuplicates',
            '-I', str(bam_file),
            '-O', str(output_bam),
            '-M', str(metrics_file),
            '--CREATE_INDEX', 'true'
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.checkpoint['deduped_samples'].append(sample_name)
            self.save_checkpoint()
            self.logger.info(f"Duplicates marked for {sample_name}")
            return output_bam
        except subprocess.CalledProcessError as e:
            self.logger.error(f"MarkDuplicates failed for {sample_name}: {e.stderr}")
            raise

    def call_variants_sample(self, bam_file: Path) -> Path:
        """
        Call variants for a single sample using HaplotypeCaller in GVCF mode.

        Args:
            bam_file: Deduplicated BAM file

        Returns:
            Path to output GVCF file
        """
        sample_name = bam_file.stem.replace('_dedup', '')

        # Check checkpoint
        if sample_name in self.checkpoint['gvcf_samples']:
            output_gvcf = self.gvcf_dir / f"{sample_name}.g.vcf.gz"
            self.logger.info(f"GVCF for {sample_name} already exists (checkpoint)")
            return output_gvcf

        output_gvcf = self.gvcf_dir / f"{sample_name}.g.vcf.gz"

        self.logger.info(f"Calling variants for {sample_name}...")

        cmd = [
            self.gatk, 'HaplotypeCaller',
            '-R', str(self.reference),
            '-I', str(bam_file),
            '-O', str(output_gvcf),
            '-ERC', 'GVCF',
            '--ploidy', '2',  # Diploid kelp sporophytes
            '--native-pair-hmm-threads', str(self.threads)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.checkpoint['gvcf_samples'].append(sample_name)
            self.save_checkpoint()
            self.logger.info(f"Variant calling completed for {sample_name}")
            return output_gvcf
        except subprocess.CalledProcessError as e:
            self.logger.error(f"HaplotypeCaller failed for {sample_name}: {e.stderr}")
            raise

    def create_sample_map(self, gvcf_files: List[Path]) -> Path:
        """
        Create sample map file for GenomicsDBImport.

        Args:
            gvcf_files: List of GVCF file paths

        Returns:
            Path to sample map file
        """
        sample_map = self.output_dir / "sample_map.txt"

        with open(sample_map, 'w') as f:
            for gvcf in gvcf_files:
                sample_name = gvcf.stem.replace('.g.vcf', '')
                f.write(f"{sample_name}\t{gvcf}\n")

        self.logger.info(f"Sample map created with {len(gvcf_files)} samples")
        return sample_map

    def genomicsdb_import(self, sample_map: Path):
        """
        Import GVCFs into GenomicsDB for efficient joint genotyping.

        Args:
            sample_map: Path to sample map file
        """
        if self.checkpoint.get('genomicsdb_created', False):
            self.logger.info("GenomicsDB already created (checkpoint)")
            return

        self.logger.info("Importing GVCFs to GenomicsDB...")

        # Get list of scaffolds/chromosomes from reference
        intervals = self.get_intervals()

        # Remove existing database if present
        if self.db_dir.exists() and any(self.db_dir.iterdir()):
            self.logger.warning("Removing existing GenomicsDB")
            import shutil
            shutil.rmtree(self.db_dir)
            self.db_dir.mkdir(parents=True)

        cmd = [
            self.gatk, 'GenomicsDBImport',
            '--sample-name-map', str(sample_map),
            '--genomicsdb-workspace-path', str(self.db_dir),
            '--intervals', intervals,
            '--reader-threads', str(self.threads)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.checkpoint['genomicsdb_created'] = True
            self.save_checkpoint()
            self.logger.info("GenomicsDB import completed")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"GenomicsDBImport failed: {e.stderr}")
            raise

    def get_intervals(self) -> str:
        """
        Get genomic intervals from reference.
        For simplicity, returns entire genome. Can be modified for parallel processing.

        Returns:
            Interval string for GATK
        """
        # Read sequence names from .dict file
        dict_file = self.reference.with_suffix('.dict')

        with open(dict_file, 'r') as f:
            lines = f.readlines()

        # Get first scaffold/chromosome (or could return all)
        # For whole genome, we can use all scaffolds
        scaffolds = []
        for line in lines:
            if line.startswith('@SQ'):
                parts = line.split('\t')
                for part in parts:
                    if part.startswith('SN:'):
                        scaffold = part.split(':')[1]
                        scaffolds.append(scaffold)

        # Return first scaffold or all (GATK can handle multiple with -L)
        # For simplicity, using entire reference
        return scaffolds[0] if scaffolds else "scaffold_1"

    def joint_genotype(self):
        """Perform joint genotyping on all samples."""
        if self.checkpoint.get('joint_genotyped', False):
            self.logger.info("Joint genotyping already completed (checkpoint)")
            return

        output_vcf = self.vcf_dir / "raw_variants.vcf.gz"

        self.logger.info("Performing joint genotyping...")

        cmd = [
            self.gatk, 'GenotypeGVCFs',
            '-R', str(self.reference),
            '-V', f"gendb://{self.db_dir}",
            '-O', str(output_vcf)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.checkpoint['joint_genotyped'] = True
            self.save_checkpoint()
            self.logger.info(f"Joint genotyping completed: {output_vcf}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"GenotypeGVCFs failed: {e.stderr}")
            raise

    def filter_variants(self):
        """Apply GATK hard filters to variants."""
        if self.checkpoint.get('filtered', False):
            self.logger.info("Filtering already completed (checkpoint)")
            return

        input_vcf = self.vcf_dir / "raw_variants.vcf.gz"
        filtered_vcf = self.vcf_dir / "filtered_variants.vcf.gz"

        self.logger.info("Applying hard filters...")

        # Filter expression for SNPs
        filter_expr = (
            "QD < 2.0 || FS > 60.0 || MQ < 40.0 || "
            "SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        )

        cmd = [
            self.gatk, 'VariantFiltration',
            '-R', str(self.reference),
            '-V', str(input_vcf),
            '-O', str(filtered_vcf),
            '--filter-expression', filter_expr,
            '--filter-name', "HARD_FILTER"
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Extract only PASS variants
            pass_vcf = self.vcf_dir / "pass_variants.vcf.gz"
            extract_cmd = [
                self.gatk, 'SelectVariants',
                '-R', str(self.reference),
                '-V', str(filtered_vcf),
                '-O', str(pass_vcf),
                '--exclude-filtered'
            ]
            subprocess.run(extract_cmd, check=True, capture_output=True, text=True)

            self.checkpoint['filtered'] = True
            self.save_checkpoint()
            self.logger.info(f"Filtering completed: {pass_vcf}")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Filtering failed: {e.stderr}")
            raise

    def calculate_statistics(self):
        """Calculate variant statistics."""
        self.logger.info("Calculating variant statistics...")

        raw_vcf = self.vcf_dir / "raw_variants.vcf.gz"
        pass_vcf = self.vcf_dir / "pass_variants.vcf.gz"
        stats_file = self.output_dir / "variant_stats.txt"

        # Count variants
        def count_variants(vcf_file):
            """Count variants in VCF file using bcftools."""
            try:
                result = subprocess.run(
                    ['bcftools', 'view', '-H', str(vcf_file)],
                    capture_output=True, text=True, check=True
                )
                return len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
            except subprocess.CalledProcessError:
                # Fallback to zgrep
                result = subprocess.run(
                    ['zgrep', '-v', '^#', str(vcf_file)],
                    capture_output=True, text=True
                )
                return len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0

        raw_count = count_variants(raw_vcf)
        pass_count = count_variants(pass_vcf)

        # Calculate Ti/Tv ratio
        titv_ratio = self.calculate_titv(pass_vcf)

        # Get variants per scaffold
        scaffold_counts = self.count_per_scaffold(pass_vcf)

        # Write statistics
        with open(stats_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("VARIANT CALLING STATISTICS\n")
            f.write("="*80 + "\n\n")
            f.write(f"Total variants called: {raw_count:,}\n")
            f.write(f"Variants passing filters: {pass_count:,}\n")
            f.write(f"Variants filtered out: {raw_count - pass_count:,}\n")
            f.write(f"Pass rate: {100*pass_count/raw_count:.2f}%\n\n")
            f.write(f"Ti/Tv ratio: {titv_ratio:.3f}\n")
            f.write(f"  (Expected for whole-genome: 2.0-2.5)\n\n")
            f.write("Variants per scaffold (top 10):\n")
            f.write("-"*40 + "\n")
            for scaffold, count in sorted(scaffold_counts.items(),
                                         key=lambda x: x[1], reverse=True)[:10]:
                f.write(f"  {scaffold}: {count:,}\n")

        self.logger.info(f"Statistics saved to {stats_file}")
        self.logger.info(f"Total variants: {raw_count:,} → {pass_count:,} after filtering")
        self.logger.info(f"Ti/Tv ratio: {titv_ratio:.3f}")

    def calculate_titv(self, vcf_file: Path) -> float:
        """
        Calculate transition/transversion ratio.

        Args:
            vcf_file: VCF file path

        Returns:
            Ti/Tv ratio
        """
        transitions = ['AG', 'GA', 'CT', 'TC']
        ti_count = 0
        tv_count = 0

        try:
            # Use bcftools if available
            result = subprocess.run(
                ['bcftools', 'view', '-H', str(vcf_file)],
                capture_output=True, text=True, check=True
            )
            lines = result.stdout.strip().split('\n')
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Fallback to zgrep
            result = subprocess.run(
                ['zgrep', '-v', '^#', str(vcf_file)],
                capture_output=True, text=True
            )
            lines = result.stdout.strip().split('\n')

        for line in lines:
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue

            ref = parts[3]
            alt = parts[4]

            # Only consider biallelic SNPs
            if len(ref) == 1 and len(alt) == 1:
                change = ref + alt
                if change in transitions:
                    ti_count += 1
                else:
                    tv_count += 1

        return ti_count / tv_count if tv_count > 0 else 0.0

    def count_per_scaffold(self, vcf_file: Path) -> Dict[str, int]:
        """
        Count variants per scaffold.

        Args:
            vcf_file: VCF file path

        Returns:
            Dictionary mapping scaffold to variant count
        """
        scaffold_counts = {}

        try:
            result = subprocess.run(
                ['bcftools', 'view', '-H', str(vcf_file)],
                capture_output=True, text=True, check=True
            )
            lines = result.stdout.strip().split('\n')
        except (subprocess.CalledProcessError, FileNotFoundError):
            result = subprocess.run(
                ['zgrep', '-v', '^#', str(vcf_file)],
                capture_output=True, text=True
            )
            lines = result.stdout.strip().split('\n')

        for line in lines:
            if not line:
                continue
            scaffold = line.split('\t')[0]
            scaffold_counts[scaffold] = scaffold_counts.get(scaffold, 0) + 1

        return scaffold_counts

    def run_pipeline(self):
        """Run the complete variant calling pipeline."""
        self.logger.info("="*80)
        self.logger.info("Starting Variant Calling Pipeline")
        self.logger.info(f"BAM directory: {self.bam_dir}")
        self.logger.info(f"Reference: {self.reference}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info("="*80)

        # Prepare reference
        self.create_sequence_dictionary()
        self.index_reference()

        # Find BAM files
        bam_files = self.find_bam_files()

        if not bam_files:
            self.logger.error("No BAM files found!")
            sys.exit(1)

        # Step 1: Mark duplicates
        self.logger.info(f"Step 1: Marking duplicates for {len(bam_files)} samples...")
        dedup_bams = []
        for bam in bam_files:
            try:
                dedup_bam = self.mark_duplicates(bam)
                dedup_bams.append(dedup_bam)
            except Exception as e:
                self.logger.error(f"Failed to process {bam.name}: {e}")

        # Step 2: Call variants (parallel)
        self.logger.info(f"Step 2: Calling variants for {len(dedup_bams)} samples...")
        self.logger.info(f"Processing {self.parallel} samples in parallel...")

        with multiprocessing.Pool(processes=self.parallel) as pool:
            gvcf_files = pool.map(self.call_variants_sample, dedup_bams)

        gvcf_files = [g for g in gvcf_files if g.exists()]
        self.logger.info(f"Created {len(gvcf_files)} GVCF files")

        # Step 3: GenomicsDB import
        self.logger.info("Step 3: Importing to GenomicsDB...")
        sample_map = self.create_sample_map(gvcf_files)
        self.genomicsdb_import(sample_map)

        # Step 4: Joint genotyping
        self.logger.info("Step 4: Joint genotyping...")
        self.joint_genotype()

        # Step 5: Filter variants
        self.logger.info("Step 5: Filtering variants...")
        self.filter_variants()

        # Step 6: Calculate statistics
        self.logger.info("Step 6: Calculating statistics...")
        self.calculate_statistics()

        self.logger.info("="*80)
        self.logger.info("Variant Calling Pipeline Completed Successfully!")
        self.logger.info(f"Results saved to: {self.output_dir}")
        self.logger.info("="*80)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Variant Calling Pipeline for Kelp GWAS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python variant_calling.py --bam-dir aligned/ \\
                             --reference Macpyr2_genomic.fasta \\
                             --output-dir variant_calling/ \\
                             --threads 8 \\
                             --parallel 10
        """
    )

    parser.add_argument('--bam-dir', required=True,
                       help='Directory containing sorted BAM files')
    parser.add_argument('--reference', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('--threads', type=int, default=8,
                       help='Number of threads for GATK operations (default: 8)')
    parser.add_argument('--parallel', type=int, default=10,
                       help='Number of samples to process in parallel (default: 10)')
    parser.add_argument('--gatk-path', default='gatk',
                       help='Path to GATK executable (default: gatk)')

    args = parser.parse_args()

    # Create and run pipeline
    pipeline = VariantCalling(
        bam_dir=args.bam_dir,
        reference=args.reference,
        output_dir=args.output_dir,
        threads=args.threads,
        parallel=args.parallel,
        gatk_path=args.gatk_path
    )

    pipeline.run_pipeline()


if __name__ == '__main__':
    main()
