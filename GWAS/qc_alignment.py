#!/usr/bin/env python3
"""
Quality Control and Alignment Module for Kelp GWAS Pipeline

This module performs:
1. FastQC quality assessment on raw FASTQ files
2. Adapter trimming and quality filtering with fastp
3. Reference genome indexing with HISAT2
4. Read alignment to reference genome
5. BAM file sorting and indexing
6. Alignment statistics and QC summary generation

Author: GWAS Pipeline
Date: 2026-01-26
Usage:
    python qc_alignment.py --fastq-dir raw_data/ --reference Macpyr2_genomic.fasta \\
                           --output-dir qc_output/ --threads 8 --parallel 4
"""

import os
import sys
import logging
import argparse
import subprocess
import multiprocessing
from pathlib import Path
from typing import List, Tuple, Dict
import pandas as pd
import shutil
from datetime import datetime
import json


class QCAlignment:
    """
    Quality control and alignment pipeline for kelp genome sequencing data.
    """

    def __init__(self, fastq_dir: str, reference: str, output_dir: str,
                 threads: int = 8, parallel: int = 4):
        """
        Initialize QC and alignment pipeline.

        Args:
            fastq_dir: Directory containing paired-end FASTQ files
            reference: Path to reference genome FASTA file
            output_dir: Output directory for results
            threads: Number of threads per alignment job
            parallel: Number of samples to process in parallel
        """
        self.fastq_dir = Path(fastq_dir)
        self.reference = Path(reference)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.parallel = parallel

        # Create output directories
        self.qc_reports_dir = self.output_dir / "qc_reports"
        self.trimmed_dir = self.output_dir / "trimmed"
        self.aligned_dir = self.output_dir / "aligned"

        for dir_path in [self.qc_reports_dir, self.trimmed_dir, self.aligned_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        # Setup logging
        self.log_file = self.output_dir / "qc_pipeline.log"
        self.setup_logging()

        # Checkpoint file
        self.checkpoint_file = self.output_dir / "qc_checkpoint.json"
        self.checkpoint = self.load_checkpoint()

        # QC metrics storage
        self.qc_metrics = []

    def setup_logging(self):
        """Configure logging with file and console handlers."""
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
        """Load checkpoint data if exists."""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        return {'completed_samples': [], 'index_built': False}

    def save_checkpoint(self):
        """Save checkpoint data."""
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.checkpoint, f, indent=2)

    def check_disk_space(self, required_gb: int = 500):
        """
        Check if sufficient disk space is available.

        Args:
            required_gb: Required space in gigabytes

        Raises:
            IOError: If insufficient disk space
        """
        stat = shutil.disk_usage(self.output_dir)
        free_gb = stat.free / (1024**3)

        if free_gb < required_gb:
            raise IOError(f"Insufficient disk space: {free_gb:.1f} GB available, "
                         f"{required_gb} GB required")

        self.logger.info(f"Disk space check: {free_gb:.1f} GB available")

    def check_dependencies(self):
        """Check if required tools are installed."""
        required_tools = ['fastqc', 'fastp', 'hisat2', 'hisat2-build', 'samtools']
        missing_tools = []

        for tool in required_tools:
            if shutil.which(tool) is None:
                missing_tools.append(tool)

        if missing_tools:
            raise EnvironmentError(f"Missing required tools: {', '.join(missing_tools)}. "
                                  f"Please install them before running the pipeline.")

        self.logger.info("All required tools are available")

    def find_fastq_pairs(self) -> List[Tuple[Path, Path, str]]:
        """
        Find paired-end FASTQ files in the input directory.

        Returns:
            List of tuples: (R1_path, R2_path, sample_name)
        """
        # Common patterns for paired-end files
        fastq_files = list(self.fastq_dir.glob("*.fastq.gz")) + \
                     list(self.fastq_dir.glob("*.fq.gz"))

        # Group by sample name
        r1_files = sorted([f for f in fastq_files if '_R1' in f.name or '_1.f' in f.name])
        pairs = []

        for r1_file in r1_files:
            # Find corresponding R2 file
            r2_file = None
            sample_name = r1_file.name.replace('_R1', '').replace('_1.f', '')

            for pattern in ['_R2', '_2.f']:
                potential_r2 = r1_file.parent / r1_file.name.replace(
                    r1_file.name.split('_')[0] + '_R1' if '_R1' in r1_file.name else r1_file.name.split('_')[0] + '_1',
                    r1_file.name.split('_')[0] + pattern
                )
                if potential_r2.exists():
                    r2_file = potential_r2
                    break

            if r2_file:
                # Extract sample name
                sample_name = r1_file.name.split('_R')[0].split('_1')[0]
                pairs.append((r1_file, r2_file, sample_name))
            else:
                self.logger.warning(f"Could not find R2 pair for {r1_file.name}")

        self.logger.info(f"Found {len(pairs)} paired-end FASTQ file sets")
        return pairs

    def run_fastqc(self, fastq_files: List[Path]) -> bool:
        """
        Run FastQC on FASTQ files.

        Args:
            fastq_files: List of FASTQ file paths

        Returns:
            True if successful
        """
        self.logger.info(f"Running FastQC on {len(fastq_files)} files...")

        cmd = [
            'fastqc',
            '-o', str(self.qc_reports_dir),
            '-t', str(self.threads)
        ] + [str(f) for f in fastq_files]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("FastQC completed successfully")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"FastQC failed: {e.stderr}")
            return False

    def run_fastp(self, r1: Path, r2: Path, sample_name: str) -> Tuple[Path, Path, Dict]:
        """
        Run fastp for adapter trimming and quality filtering.

        Args:
            r1: Read 1 FASTQ file
            r2: Read 2 FASTQ file
            sample_name: Sample identifier

        Returns:
            Tuple of (trimmed_r1, trimmed_r2, metrics_dict)
        """
        out_r1 = self.trimmed_dir / f"{sample_name}_R1_trimmed.fastq.gz"
        out_r2 = self.trimmed_dir / f"{sample_name}_R2_trimmed.fastq.gz"
        json_report = self.trimmed_dir / f"{sample_name}_fastp.json"
        html_report = self.trimmed_dir / f"{sample_name}_fastp.html"

        cmd = [
            'fastp',
            '-i', str(r1),
            '-I', str(r2),
            '-o', str(out_r1),
            '-O', str(out_r2),
            '--detect_adapter_for_pe',
            '--qualified_quality_phred', '20',
            '--length_required', '50',
            '--thread', str(self.threads),
            '--json', str(json_report),
            '--html', str(html_report)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Parse metrics from JSON
            with open(json_report, 'r') as f:
                report = json.load(f)

            metrics = {
                'sample': sample_name,
                'total_reads_before': report['summary']['before_filtering']['total_reads'],
                'total_reads_after': report['summary']['after_filtering']['total_reads'],
                'q30_rate_before': report['summary']['before_filtering'].get('q30_rate', 0),
                'q30_rate_after': report['summary']['after_filtering'].get('q30_rate', 0),
            }

            self.logger.info(f"fastp completed for {sample_name}: "
                           f"{metrics['total_reads_after']:,} reads retained")

            return out_r1, out_r2, metrics

        except subprocess.CalledProcessError as e:
            self.logger.error(f"fastp failed for {sample_name}: {e.stderr}")
            raise

    def build_hisat2_index(self) -> bool:
        """
        Build HISAT2 index from reference genome.

        Returns:
            True if successful or already exists
        """
        index_base = self.output_dir / "reference_index" / "genome"
        index_base.parent.mkdir(exist_ok=True)

        # Check if index already exists
        if (index_base.parent / "genome.1.ht2").exists():
            self.logger.info("HISAT2 index already exists, skipping build")
            return True

        if self.checkpoint.get('index_built', False):
            self.logger.info("Index already built (from checkpoint)")
            return True

        self.logger.info(f"Building HISAT2 index from {self.reference}...")

        cmd = [
            'hisat2-build',
            '-p', str(self.threads),
            str(self.reference),
            str(index_base)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("HISAT2 index built successfully")
            self.checkpoint['index_built'] = True
            self.save_checkpoint()
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"HISAT2 index build failed: {e.stderr}")
            return False

    def align_sample(self, r1: Path, r2: Path, sample_name: str) -> Tuple[Path, Dict]:
        """
        Align reads to reference genome using HISAT2.

        Args:
            r1: Trimmed read 1 file
            r2: Trimmed read 2 file
            sample_name: Sample identifier

        Returns:
            Tuple of (sorted_bam_path, alignment_metrics)
        """
        index_base = self.output_dir / "reference_index" / "genome"
        sam_file = self.aligned_dir / f"{sample_name}.sam"
        bam_file = self.aligned_dir / f"{sample_name}.bam"
        sorted_bam = self.aligned_dir / f"{sample_name}_sorted.bam"

        # HISAT2 alignment
        self.logger.info(f"Aligning {sample_name}...")
        hisat2_cmd = [
            'hisat2',
            '-x', str(index_base),
            '-1', str(r1),
            '-2', str(r2),
            '-S', str(sam_file),
            '-p', str(self.threads),
            '--summary-file', str(self.aligned_dir / f"{sample_name}_alignment_summary.txt")
        ]

        try:
            result = subprocess.run(hisat2_cmd, check=True, capture_output=True, text=True)

            # Parse alignment rate from stderr
            alignment_rate = 0.0
            for line in result.stderr.split('\n'):
                if 'overall alignment rate' in line:
                    alignment_rate = float(line.split('%')[0].strip())

            # Convert SAM to BAM
            subprocess.run(['samtools', 'view', '-bS', str(sam_file), '-o', str(bam_file)],
                         check=True, capture_output=True)

            # Sort BAM
            subprocess.run(['samtools', 'sort', '-@', str(self.threads),
                          str(bam_file), '-o', str(sorted_bam)],
                         check=True, capture_output=True)

            # Index BAM
            subprocess.run(['samtools', 'index', str(sorted_bam)],
                         check=True, capture_output=True)

            # Calculate coverage statistics
            stats_result = subprocess.run(
                ['samtools', 'coverage', str(sorted_bam)],
                check=True, capture_output=True, text=True
            )

            # Parse coverage
            mean_depth = 0.0
            breadth_4x = 0.0
            for line in stats_result.stdout.split('\n')[1:]:  # Skip header
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 7:
                        mean_depth += float(parts[6])  # meandepth column

            # Get detailed stats
            flagstat_result = subprocess.run(
                ['samtools', 'flagstat', str(sorted_bam)],
                check=True, capture_output=True, text=True
            )

            # Calculate breadth of coverage at 4x
            depth_result = subprocess.run(
                ['samtools', 'depth', str(sorted_bam)],
                capture_output=True, text=True
            )

            total_positions = 0
            covered_4x = 0
            for line in depth_result.stdout.split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        total_positions += 1
                        if int(parts[2]) >= 4:
                            covered_4x += 1

            if total_positions > 0:
                breadth_4x = (covered_4x / total_positions) * 100

            # Clean up intermediate files
            sam_file.unlink()
            bam_file.unlink()

            metrics = {
                'sample': sample_name,
                'alignment_rate': alignment_rate,
                'mean_depth': mean_depth,
                'breadth_4x': breadth_4x
            }

            # Check QC thresholds
            qc_pass = alignment_rate >= 80.0 and mean_depth >= 4.0
            metrics['qc_pass'] = qc_pass

            if not qc_pass:
                self.logger.warning(f"Sample {sample_name} failed QC: "
                                  f"alignment_rate={alignment_rate:.1f}%, "
                                  f"mean_depth={mean_depth:.1f}x")

            self.logger.info(f"Alignment completed for {sample_name}: "
                           f"{alignment_rate:.1f}% aligned, "
                           f"{mean_depth:.1f}x mean depth")

            return sorted_bam, metrics

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Alignment failed for {sample_name}: {e.stderr}")
            raise

    def process_sample(self, r1: Path, r2: Path, sample_name: str) -> Dict:
        """
        Process a single sample through the QC and alignment pipeline.

        Args:
            r1: Read 1 FASTQ file
            r2: Read 2 FASTQ file
            sample_name: Sample identifier

        Returns:
            Dictionary with QC metrics
        """
        # Check if already completed
        if sample_name in self.checkpoint['completed_samples']:
            self.logger.info(f"Sample {sample_name} already processed (from checkpoint)")
            return {}

        try:
            # Run fastp
            trimmed_r1, trimmed_r2, trim_metrics = self.run_fastp(r1, r2, sample_name)

            # Run alignment
            sorted_bam, align_metrics = self.align_sample(trimmed_r1, trimmed_r2, sample_name)

            # Combine metrics
            metrics = {**trim_metrics, **align_metrics}

            # Update checkpoint
            self.checkpoint['completed_samples'].append(sample_name)
            self.save_checkpoint()

            return metrics

        except Exception as e:
            self.logger.error(f"Failed to process {sample_name}: {str(e)}")
            return {'sample': sample_name, 'error': str(e), 'qc_pass': False}

    def run_pipeline(self):
        """
        Run the complete QC and alignment pipeline.
        """
        self.logger.info("="*80)
        self.logger.info("Starting QC and Alignment Pipeline")
        self.logger.info(f"Input directory: {self.fastq_dir}")
        self.logger.info(f"Reference: {self.reference}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Threads per job: {self.threads}")
        self.logger.info(f"Parallel jobs: {self.parallel}")
        self.logger.info("="*80)

        # Check dependencies and disk space
        self.check_dependencies()
        self.check_disk_space()

        # Find FASTQ pairs
        fastq_pairs = self.find_fastq_pairs()

        if not fastq_pairs:
            self.logger.error("No paired-end FASTQ files found!")
            sys.exit(1)

        # Run FastQC on all raw files (only once)
        all_fastq_files = []
        for r1, r2, _ in fastq_pairs:
            all_fastq_files.extend([r1, r2])
        self.run_fastqc(all_fastq_files)

        # Build HISAT2 index
        if not self.build_hisat2_index():
            self.logger.error("Failed to build HISAT2 index")
            sys.exit(1)

        # Process samples in parallel
        self.logger.info(f"Processing {len(fastq_pairs)} samples with {self.parallel} parallel jobs...")

        # Use multiprocessing pool
        with multiprocessing.Pool(processes=self.parallel) as pool:
            results = pool.starmap(self.process_sample, fastq_pairs)

        # Collect all metrics
        self.qc_metrics = [r for r in results if r]  # Filter empty results

        # Generate summary
        self.generate_summary()

        self.logger.info("="*80)
        self.logger.info("QC and Alignment Pipeline Completed Successfully!")
        self.logger.info(f"Results saved to: {self.output_dir}")
        self.logger.info("="*80)

    def generate_summary(self):
        """Generate QC summary report."""
        if not self.qc_metrics:
            self.logger.warning("No QC metrics to summarize")
            return

        # Create DataFrame
        df = pd.DataFrame(self.qc_metrics)

        # Save detailed metrics
        summary_file = self.output_dir / "qc_summary.csv"
        df.to_csv(summary_file, index=False)

        # Generate summary statistics
        summary_text = self.output_dir / "qc_summary.txt"
        with open(summary_text, 'w') as f:
            f.write("="*80 + "\n")
            f.write("QC AND ALIGNMENT SUMMARY\n")
            f.write("="*80 + "\n\n")

            f.write(f"Total samples processed: {len(df)}\n")
            f.write(f"Samples passing QC: {df['qc_pass'].sum()}\n")
            f.write(f"Samples failing QC: {(~df['qc_pass']).sum()}\n\n")

            if 'alignment_rate' in df.columns:
                f.write(f"Mean alignment rate: {df['alignment_rate'].mean():.2f}%\n")
                f.write(f"Mean coverage depth: {df['mean_depth'].mean():.2f}x\n")
                f.write(f"Mean breadth at 4x: {df['breadth_4x'].mean():.2f}%\n\n")

            # List failed samples
            failed = df[~df['qc_pass']]
            if len(failed) > 0:
                f.write("FAILED SAMPLES:\n")
                f.write("-"*80 + "\n")
                for _, row in failed.iterrows():
                    f.write(f"  {row['sample']}: ")
                    if 'error' in row:
                        f.write(f"Error - {row['error']}\n")
                    else:
                        f.write(f"alignment={row.get('alignment_rate', 0):.1f}%, "
                               f"depth={row.get('mean_depth', 0):.1f}x\n")

        self.logger.info(f"QC summary saved to {summary_file}")


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Quality Control and Alignment Pipeline for Kelp GWAS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python qc_alignment.py --fastq-dir raw_data/ \\
                          --reference Macpyr2_genomic.fasta \\
                          --output-dir qc_output/ \\
                          --threads 8 \\
                          --parallel 4
        """
    )

    parser.add_argument('--fastq-dir', required=True,
                       help='Directory containing paired-end FASTQ files')
    parser.add_argument('--reference', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for results')
    parser.add_argument('--threads', type=int, default=8,
                       help='Number of threads per alignment job (default: 8)')
    parser.add_argument('--parallel', type=int, default=4,
                       help='Number of samples to process in parallel (default: 4)')

    args = parser.parse_args()

    # Create and run pipeline
    pipeline = QCAlignment(
        fastq_dir=args.fastq_dir,
        reference=args.reference,
        output_dir=args.output_dir,
        threads=args.threads,
        parallel=args.parallel
    )

    pipeline.run_pipeline()


if __name__ == '__main__':
    main()
