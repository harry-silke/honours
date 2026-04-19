#!/usr/bin/env python3
"""
Main Pipeline Orchestrator for Kelp GWAS Analysis

This script coordinates the execution of all pipeline modules:
1. Quality control and alignment
2. Variant calling
3. VCF filtering
4. Population structure analysis (R)
5. GWAS analysis (R)
6. Visualization (R)
7. Gene annotation
8. Genomic prediction (R)
9. Final report generation

Author: GWAS Pipeline
Date: 2026-01-26
Usage:
    python run_gwas_pipeline.py --config config.yaml
    python run_gwas_pipeline.py --config config.yaml --start-from filtering
"""

import os
import sys
import logging
import argparse
import subprocess
import yaml
import json
from pathlib import Path
from datetime import datetime
import shutil


class GWASPipeline:
    """Main GWAS pipeline orchestrator."""

    def __init__(self, config_file: str, start_from: str = None):
        """
        Initialize pipeline.

        Args:
            config_file: Path to YAML configuration file
            start_from: Pipeline step to start from (for resuming)
        """
        self.config_file = Path(config_file)
        self.start_from = start_from

        # Load configuration
        self.load_config()

        # Setup output directory
        self.base_output = Path(self.config['output']['base_dir'])
        self.base_output.mkdir(parents=True, exist_ok=True)

        # Setup logging
        self.setup_logging()

        # Checkpoint file
        self.checkpoint_file = self.base_output / "pipeline_checkpoint.json"
        self.checkpoint = self.load_checkpoint()

        # Define pipeline steps
        self.steps = [
            'qc',
            'variant_calling',
            'filtering',
            'population_structure',
            'gwas',
            'visualization',
            'annotation',
            'genomic_prediction'
        ]

    def load_config(self):
        """Load configuration from YAML file."""
        with open(self.config_file, 'r') as f:
            self.config = yaml.safe_load(f)

    def setup_logging(self):
        """Configure logging."""
        log_file = self.base_output / "pipeline.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def load_checkpoint(self):
        """Load checkpoint data."""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        return {'completed_steps': [], 'failed_steps': []}

    def save_checkpoint(self):
        """Save checkpoint data."""
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.checkpoint, f, indent=2)

    def run_command(self, cmd: list, description: str) -> bool:
        """
        Run a shell command.

        Args:
            cmd: Command as list of strings
            description: Description of the command

        Returns:
            True if successful, False otherwise
        """
        self.logger.info(f"Running: {description}")
        self.logger.info(f"Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            if result.stdout:
                self.logger.debug(result.stdout)

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {description}")
            self.logger.error(f"Error: {e.stderr}")
            return False

    def step_qc_alignment(self):
        """Step 1: Quality control and alignment."""
        step_name = 'qc'

        if not self.config['pipeline'].get('run_qc', True):
            self.logger.info("Skipping QC step (disabled in config)")
            return True

        output_dir = self.base_output / self.config['output']['qc']

        cmd = [
            'python', 'qc_alignment.py',
            '--fastq-dir', self.config['input']['fastq_dir'],
            '--reference', self.config['input']['reference_genome'],
            '--output-dir', str(output_dir),
            '--threads', str(self.config['resources']['threads']),
            '--parallel', str(self.config['resources']['parallel_samples'])
        ]

        return self.run_command(cmd, "Quality Control and Alignment")

    def step_variant_calling(self):
        """Step 2: Variant calling."""
        step_name = 'variant_calling'

        if not self.config['pipeline'].get('run_variant_calling', True):
            self.logger.info("Skipping variant calling (disabled in config)")
            return True

        bam_dir = self.base_output / self.config['output']['qc'] / 'aligned'
        output_dir = self.base_output / self.config['output']['variants']

        cmd = [
            'python', 'variant_calling.py',
            '--bam-dir', str(bam_dir),
            '--reference', self.config['input']['reference_genome'],
            '--output-dir', str(output_dir),
            '--threads', str(self.config['resources']['threads']),
            '--parallel', str(self.config['resources']['parallel_samples']),
            '--gatk-path', self.config['variant_calling'].get('gatk_path', 'gatk')
        ]

        return self.run_command(cmd, "Variant Calling")

    def step_filtering(self):
        """Step 3: VCF filtering."""
        step_name = 'filtering'

        if not self.config['pipeline'].get('run_filtering', True):
            self.logger.info("Skipping filtering (disabled in config)")
            return True

        vcf_file = self.base_output / self.config['output']['variants'] / 'vcfs' / 'pass_variants.vcf.gz'
        output_dir = self.base_output / self.config['output']['filtering']

        cmd = [
            'python', 'filter_vcf.py',
            '--vcf', str(vcf_file),
            '--reference', self.config['input']['reference_genome'],
            '--annotation', self.config['input']['annotation'],
            '--output-dir', str(output_dir),
            '--maf', str(self.config['vcf_filtering']['min_maf']),
            '--missing', str(self.config['vcf_filtering']['max_missing']),
            '--hwe-pval', str(self.config['vcf_filtering']['hwe_pval']),
            '--min-depth', str(self.config['vcf_filtering']['min_snp_depth']),
            '--max-depth', str(self.config['vcf_filtering']['max_snp_depth'])
        ]

        return self.run_command(cmd, "VCF Filtering")

    def step_population_structure(self):
        """Step 4: Population structure analysis (R)."""
        step_name = 'population_structure'

        if not self.config['pipeline'].get('run_population_structure', True):
            self.logger.info("Skipping population structure (disabled in config)")
            return True

        plink_prefix = self.base_output / self.config['output']['filtering'] / 'plink_format'
        output_dir = self.base_output / self.config['output']['population']

        cmd = [
            'Rscript', 'population_structure.R',
            '--plink', str(plink_prefix),
            '--phenotypes', self.config['input']['phenotypes'],
            '--output-dir', str(output_dir),
            '--n-pcs', str(self.config['population_structure']['n_pcs'])
        ]

        return self.run_command(cmd, "Population Structure Analysis")

    def step_gwas(self):
        """Step 5: GWAS analysis (R)."""
        step_name = 'gwas'

        if not self.config['pipeline'].get('run_gwas', True):
            self.logger.info("Skipping GWAS (disabled in config)")
            return True

        plink_prefix = self.base_output / self.config['output']['filtering'] / 'plink_format'
        pcs_file = self.base_output / self.config['output']['population'] / 'pca_results.csv'
        kinship_file = self.base_output / self.config['output']['population'] / 'kinship_matrix.txt'
        output_dir = self.base_output / self.config['output']['gwas']

        cmd = [
            'Rscript', 'gwas_analysis.R',
            '--geno', str(plink_prefix),
            '--pheno', self.config['input']['phenotypes'],
            '--pcs', str(pcs_file),
            '--kinship', str(kinship_file),
            '--output-dir', str(output_dir),
            '--n-pcs', str(self.config['gwas']['n_pcs'])
        ]

        return self.run_command(cmd, "GWAS Analysis")

    def step_visualization(self):
        """Step 6: GWAS visualization (R)."""
        step_name = 'visualization'

        if not self.config['pipeline'].get('run_visualization', True):
            self.logger.info("Skipping visualization (disabled in config)")
            return True

        gwas_file = self.base_output / self.config['output']['gwas'] / 'gwas_all_snps.txt'
        output_dir = self.base_output / self.config['output']['figures']

        cmd = [
            'Rscript', 'visualize_gwas.R',
            '--gwas', str(gwas_file),
            '--pheno', self.config['input']['phenotypes'],
            '--output', str(output_dir)
        ]

        return self.run_command(cmd, "GWAS Visualization")

    def step_annotation(self):
        """Step 7: Gene annotation."""
        step_name = 'annotation'

        if not self.config['pipeline'].get('run_annotation', True):
            self.logger.info("Skipping annotation (disabled in config)")
            return True

        # Check if significant SNPs exist
        sig_snps = self.base_output / self.config['output']['gwas'] / 'gwas_significant.txt'

        if not sig_snps.exists():
            self.logger.warning("No significant SNPs found, using top 100 instead")
            sig_snps = self.base_output / self.config['output']['gwas'] / 'gwas_top100.txt'

        output_dir = self.base_output / self.config['output']['annotation']

        cmd = [
            'python', 'annotate_genes.py',
            '--snps', str(sig_snps),
            '--annotation', self.config['input']['annotation'],
            '--output-dir', str(output_dir),
            '--window', str(self.config['annotation']['gene_window'])
        ]

        return self.run_command(cmd, "Gene Annotation")

    def step_genomic_prediction(self):
        """Step 8: Genomic prediction (R)."""
        step_name = 'genomic_prediction'

        if not self.config['pipeline'].get('run_genomic_prediction', True):
            self.logger.info("Skipping genomic prediction (disabled in config)")
            return True

        plink_prefix = self.base_output / self.config['output']['filtering'] / 'plink_format'
        gwas_file = self.base_output / self.config['output']['gwas'] / 'gwas_all_snps.txt'
        output_dir = self.base_output / self.config['output']['prediction']

        cmd = [
            'Rscript', 'genomic_prediction.R',
            '--geno', str(plink_prefix),
            '--pheno', self.config['input']['phenotypes'],
            '--gwas', str(gwas_file),
            '--output', str(output_dir)
        ]

        return self.run_command(cmd, "Genomic Prediction")

    def generate_final_report(self):
        """Generate comprehensive final report."""
        self.logger.info("Generating final report...")

        report_file = self.base_output / "final_report.txt"

        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("KELP GWAS PIPELINE - FINAL REPORT\n")
            f.write("="*80 + "\n\n")

            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Config File: {self.config_file}\n\n")

            # Completed steps
            f.write("Pipeline Steps:\n")
            f.write("-"*80 + "\n")
            for step in self.steps:
                status = "✓ COMPLETED" if step in self.checkpoint['completed_steps'] else "✗ NOT RUN"
                f.write(f"  {step:30s} {status}\n")
            f.write("\n")

            # Key results
            f.write("Key Results:\n")
            f.write("-"*80 + "\n")

            # Try to read key metrics
            try:
                # Number of samples
                qc_summary = self.base_output / self.config['output']['qc'] / 'qc_summary.csv'
                if qc_summary.exists():
                    import pandas as pd
                    qc_df = pd.read_csv(qc_summary)
                    f.write(f"  Samples processed: {len(qc_df)}\n")

                # Number of SNPs
                snp_stats = self.base_output / self.config['output']['filtering'] / 'snp_statistics.csv'
                if snp_stats.exists():
                    snp_df = pd.read_csv(snp_stats)
                    f.write(f"  High-quality SNPs: {len(snp_df):,}\n")

                # Significant SNPs
                sig_snps = self.base_output / self.config['output']['gwas'] / 'gwas_significant.txt'
                if sig_snps.exists():
                    sig_df = pd.read_csv(sig_snps, sep='\t')
                    f.write(f"  Significant SNPs: {len(sig_df)}\n")

                # Genomic inflation
                lambda_file = self.base_output / self.config['output']['gwas'] / 'lambda.txt'
                if lambda_file.exists():
                    with open(lambda_file) as lf:
                        lambda_val = float(lf.read().strip())
                        f.write(f"  Genomic inflation (λ): {lambda_val:.3f}\n")

                # Prediction accuracy
                model_comp = self.base_output / self.config['output']['prediction'] / 'model_comparison.csv'
                if model_comp.exists():
                    model_df = pd.read_csv(model_comp)
                    best_auc = model_df['mean_auc'].max()
                    f.write(f"  Best prediction AUC: {best_auc:.3f}\n")

            except Exception as e:
                f.write(f"  (Could not extract all metrics: {e})\n")

            f.write("\n")

            # Output files
            f.write("Output Directory Structure:\n")
            f.write("-"*80 + "\n")
            f.write(f"  {self.base_output}/\n")
            for step_dir in self.base_output.iterdir():
                if step_dir.is_dir():
                    f.write(f"    ├── {step_dir.name}/\n")
            f.write("\n")

            # Recommendations
            f.write("Next Steps:\n")
            f.write("-"*80 + "\n")
            f.write("  1. Review Manhattan plot: " +
                   str(self.base_output / self.config['output']['figures'] / 'manhattan_plot.pdf') + "\n")
            f.write("  2. Check significant SNPs: " +
                   str(self.base_output / self.config['output']['gwas'] / 'gwas_significant.txt') + "\n")
            f.write("  3. Review candidate genes: " +
                   str(self.base_output / self.config['output']['annotation'] / 'priority_genes.csv') + "\n")
            f.write("  4. Examine breeding values: " +
                   str(self.base_output / self.config['output']['prediction'] / 'gebvs_all_kelps.csv') + "\n")

            f.write("\n")
            f.write("="*80 + "\n")
            f.write("Analysis complete!\n")
            f.write("="*80 + "\n")

        self.logger.info(f"Final report saved to: {report_file}")

    def run_pipeline(self):
        """Execute the complete pipeline."""
        self.logger.info("="*80)
        self.logger.info("KELP GWAS PIPELINE")
        self.logger.info("="*80)
        self.logger.info(f"Configuration: {self.config_file}")
        self.logger.info(f"Output directory: {self.base_output}")
        self.logger.info(f"Start time: {datetime.now()}")
        self.logger.info("="*80)

        # Save copy of config
        config_copy = self.base_output / "config_used.yaml"
        shutil.copy(self.config_file, config_copy)
        self.logger.info(f"Configuration saved to: {config_copy}")

        # Determine starting point
        if self.start_from:
            if self.start_from not in self.steps:
                self.logger.error(f"Invalid start point: {self.start_from}")
                self.logger.error(f"Valid options: {', '.join(self.steps)}")
                sys.exit(1)

            start_idx = self.steps.index(self.start_from)
            steps_to_run = self.steps[start_idx:]
            self.logger.info(f"Resuming from: {self.start_from}")
        else:
            steps_to_run = self.steps

        # Define step functions
        step_functions = {
            'qc': self.step_qc_alignment,
            'variant_calling': self.step_variant_calling,
            'filtering': self.step_filtering,
            'population_structure': self.step_population_structure,
            'gwas': self.step_gwas,
            'visualization': self.step_visualization,
            'annotation': self.step_annotation,
            'genomic_prediction': self.step_genomic_prediction
        }

        # Run pipeline steps
        for step in steps_to_run:
            # Skip if already completed (from checkpoint)
            if step in self.checkpoint['completed_steps']:
                self.logger.info(f"Step '{step}' already completed (checkpoint), skipping")
                continue

            self.logger.info(f"\n{'='*80}")
            self.logger.info(f"STEP: {step.upper()}")
            self.logger.info(f"{'='*80}\n")

            # Run step
            success = step_functions[step]()

            if success:
                self.logger.info(f"✓ Step '{step}' completed successfully\n")
                self.checkpoint['completed_steps'].append(step)
                self.save_checkpoint()
            else:
                self.logger.error(f"✗ Step '{step}' FAILED\n")
                self.checkpoint['failed_steps'].append(step)
                self.save_checkpoint()

                self.logger.error("Pipeline stopped due to error")
                self.logger.error(f"To resume, run: python {sys.argv[0]} --config {self.config_file} --start-from {step}")
                sys.exit(1)

        # Generate final report
        if self.config['pipeline'].get('generate_report', True):
            self.generate_final_report()

        # Final summary
        self.logger.info("\n" + "="*80)
        self.logger.info("PIPELINE COMPLETED SUCCESSFULLY!")
        self.logger.info("="*80)
        self.logger.info(f"End time: {datetime.now()}")
        self.logger.info(f"Results saved to: {self.base_output}")
        self.logger.info(f"Final report: {self.base_output / 'final_report.txt'}")
        self.logger.info("="*80 + "\n")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Kelp GWAS Pipeline - Main Orchestrator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    # Run complete pipeline
    python run_gwas_pipeline.py --config config.yaml

    # Resume from specific step
    python run_gwas_pipeline.py --config config.yaml --start-from filtering

    # Available steps:
    qc, variant_calling, filtering, population_structure,
    gwas, visualization, annotation, genomic_prediction
        """
    )

    parser.add_argument('--config', required=True,
                       help='Path to configuration YAML file')
    parser.add_argument('--start-from', default=None,
                       help='Pipeline step to start from (for resuming)')

    args = parser.parse_args()

    # Check config file exists
    if not Path(args.config).exists():
        print(f"Error: Configuration file not found: {args.config}")
        sys.exit(1)

    # Create and run pipeline
    pipeline = GWASPipeline(
        config_file=args.config,
        start_from=args.start_from
    )

    pipeline.run_pipeline()


if __name__ == '__main__':
    main()
