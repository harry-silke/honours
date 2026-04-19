#!/usr/bin/env python3
"""
Candidate Gene Annotation Module for Kelp GWAS Pipeline

This module annotates significant SNPs with biological information:
1. Find genes near significant SNPs
2. Categorize SNPs by genomic location
3. Extract gene functions and descriptions
4. Test for functional enrichment
5. Prioritize candidate genes

Author: GWAS Pipeline
Date: 2026-01-26
Usage:
    python annotate_genes.py --snps gwas_significant.txt \\
                            --annotation Macpyr2_annotation.gff \\
                            --output-dir annotation/
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


class GeneAnnotator:
    """
    Annotate GWAS SNPs with candidate gene information.
    """

    def __init__(self, snps_file: str, annotation_file: str, output_dir: str,
                 window_size: int = 50000):
        """
        Initialize gene annotator.

        Args:
            snps_file: File with significant SNPs (TSV/CSV with chr, pos, p-value)
            annotation_file: GFF annotation file
            output_dir: Output directory
            window_size: Window size for finding nearby genes (bp)
        """
        self.snps_file = Path(snps_file)
        self.annotation_file = Path(annotation_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.window_size = window_size

        # Setup logging
        self.log_file = self.output_dir / "annotation.log"
        self.setup_logging()

        # Data storage
        self.snps = []
        self.genes = []
        self.snp_annotations = []
        self.functional_categories = self.define_functional_categories()

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

    def define_functional_categories(self) -> Dict[str, List[str]]:
        """
        Define functional gene categories and their keywords.

        Returns:
            Dictionary mapping category names to keyword lists
        """
        return {
            'Heat Shock & Stress Response': [
                'hsp', 'heat shock', 'chaperone', 'dnaj', 'dnk',
                'stress', 'groel', 'groesh', 'trigger factor'
            ],
            'Photosynthesis': [
                'photosystem', 'psii', 'psi', 'psb', 'psa', 'chlorophyll',
                'photosynthetic', 'ribulose', 'rubisco', 'calvin',
                'thylakoid', 'light-harvesting', 'lhc'
            ],
            'Antioxidant Defense': [
                'superoxide', 'dismutase', 'sod', 'catalase', 'cat',
                'peroxidase', 'pod', 'glutathione', 'ascorbate',
                'thioredoxin', 'reactive oxygen', 'ros'
            ],
            'Protein Degradation': [
                'ubiquitin', 'proteasome', 'protease', 'peptidase',
                'degradation', 'proteolysis', 'clp'
            ],
            'Membrane & Lipid': [
                'desaturase', 'lipid', 'membrane', 'fatty acid',
                'phospholipid', 'glycerol', 'acyltransferase'
            ],
            'Signal Transduction': [
                'kinase', 'receptor', 'signal', 'transduction',
                'phosphorylation', 'transcription factor'
            ],
            'Cell Wall & Structure': [
                'cell wall', 'cellulose', 'pectin', 'mannitol',
                'alginate', 'structural', 'cytoskeleton'
            ]
        }

    def load_snps(self):
        """Load significant SNPs from file."""
        self.logger.info(f"Loading SNPs from {self.snps_file}...")

        # Try to determine file format
        if self.snps_file.suffix == '.csv':
            df = pd.read_csv(self.snps_file)
        else:
            df = pd.read_csv(self.snps_file, sep='\t')

        # Normalize column names
        df.columns = df.columns.str.lower()

        # Map common column name variations
        column_map = {
            'chromosome': 'chr',
            'scaffold': 'chr',
            'chrom': 'chr',
            'pos': 'position',
            'bp': 'position',
            'p_value': 'p',
            'pvalue': 'p',
            'p.value': 'p'
        }

        for old_name, new_name in column_map.items():
            if old_name in df.columns and new_name not in df.columns:
                df.rename(columns={old_name: new_name}, inplace=True)

        # Ensure required columns exist
        required = ['chr', 'position']
        missing = [col for col in required if col not in df.columns]

        if missing:
            raise ValueError(f"SNP file missing required columns: {missing}")

        self.snps = df.to_dict('records')
        self.logger.info(f"Loaded {len(self.snps)} SNPs")

    def parse_gff_attributes(self, attr_string: str) -> Dict[str, str]:
        """
        Parse GFF attributes field.

        Args:
            attr_string: GFF attributes string

        Returns:
            Dictionary of attributes
        """
        attributes = {}
        for item in attr_string.split(';'):
            if '=' in item:
                key, value = item.strip().split('=', 1)
                attributes[key] = value
        return attributes

    def load_annotation(self):
        """Load gene annotations from GFF file."""
        self.logger.info(f"Loading annotations from {self.annotation_file}...")

        genes_loaded = 0

        with open(self.annotation_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]

                # Only process gene features
                if feature_type not in ['gene', 'mRNA', 'CDS']:
                    continue

                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = self.parse_gff_attributes(parts[8])

                # Extract gene information
                gene_id = attributes.get('ID', attributes.get('gene_id', ''))
                gene_name = attributes.get('Name', attributes.get('gene_name', gene_id))
                description = attributes.get('description', attributes.get('product', ''))
                go_terms = attributes.get('Ontology_term', '')

                if gene_id:
                    self.genes.append({
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'description': description,
                        'go_terms': go_terms,
                        'chr': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'feature_type': feature_type
                    })
                    genes_loaded += 1

        self.logger.info(f"Loaded {genes_loaded} gene features")

        # Create DataFrame for easier processing
        self.genes_df = pd.DataFrame(self.genes)

    def categorize_snp_location(self, snp_pos: int, gene_start: int, gene_end: int,
                                 strand: str, distance: int) -> str:
        """
        Categorize SNP location relative to gene.

        Args:
            snp_pos: SNP position
            gene_start: Gene start position
            gene_end: Gene end position
            strand: Gene strand (+/-)
            distance: Distance to gene (0 if within gene)

        Returns:
            Location category string
        """
        if distance == 0:
            # Within gene
            return 'Intronic/Exonic'
        elif distance <= 2000:
            # Check if upstream or downstream
            if strand == '+':
                if snp_pos < gene_start:
                    return 'Promoter'
                else:
                    return "3' UTR"
            else:
                if snp_pos > gene_end:
                    return 'Promoter'
                else:
                    return "3' UTR"
        elif distance <= self.window_size:
            return 'Regulatory Region'
        else:
            return 'Intergenic'

    def annotate_snps(self):
        """Annotate each SNP with nearest gene information."""
        self.logger.info("Annotating SNPs with gene information...")

        for snp in self.snps:
            snp_chr = str(snp['chr'])
            snp_pos = int(snp['position'])

            # Find genes on same scaffold
            same_chr_genes = self.genes_df[self.genes_df['chr'] == snp_chr]

            if len(same_chr_genes) == 0:
                # No genes on this scaffold
                self.snp_annotations.append({
                    'chr': snp_chr,
                    'position': snp_pos,
                    'gene_id': 'intergenic',
                    'gene_name': 'intergenic',
                    'distance': -1,
                    'location_type': 'Intergenic',
                    'description': '',
                    'functional_category': 'Unknown',
                    **snp
                })
                continue

            # Find nearest gene
            min_distance = float('inf')
            nearest_gene = None

            for _, gene in same_chr_genes.iterrows():
                # Check if SNP is within gene
                if gene['start'] <= snp_pos <= gene['end']:
                    distance = 0
                elif snp_pos < gene['start']:
                    distance = gene['start'] - snp_pos
                else:
                    distance = snp_pos - gene['end']

                if distance < min_distance:
                    min_distance = distance
                    nearest_gene = gene

            if nearest_gene is not None:
                # Categorize location
                location = self.categorize_snp_location(
                    snp_pos, nearest_gene['start'], nearest_gene['end'],
                    nearest_gene['strand'], min_distance
                )

                # Determine functional category
                func_category = self.categorize_gene_function(nearest_gene['description'])

                # Calculate effect size if available
                effect_size = snp.get('or', snp.get('odds_ratio', snp.get('beta', 1.0)))

                self.snp_annotations.append({
                    'chr': snp_chr,
                    'position': snp_pos,
                    'gene_id': nearest_gene['gene_id'],
                    'gene_name': nearest_gene['gene_name'],
                    'distance': min_distance,
                    'location_type': location,
                    'description': nearest_gene['description'],
                    'go_terms': nearest_gene.get('go_terms', ''),
                    'functional_category': func_category,
                    'effect_size': effect_size,
                    **snp
                })

        self.logger.info(f"Annotated {len(self.snp_annotations)} SNPs")

    def categorize_gene_function(self, description: str) -> str:
        """
        Categorize gene by function based on description.

        Args:
            description: Gene description string

        Returns:
            Functional category
        """
        description_lower = description.lower()

        for category, keywords in self.functional_categories.items():
            for keyword in keywords:
                if keyword in description_lower:
                    return category

        return 'Unknown/Other'

    def test_enrichment(self):
        """Test for functional enrichment using Fisher's exact test."""
        self.logger.info("Testing for functional enrichment...")

        annotations_df = pd.DataFrame(self.snp_annotations)

        # Count significant genes by category
        sig_counts = annotations_df['functional_category'].value_counts()

        # Count all genes by category
        all_counts = self.genes_df.apply(
            lambda x: self.categorize_gene_function(x['description']), axis=1
        ).value_counts()

        # Perform Fisher's exact test for each category
        enrichment_results = []
        total_sig = len(annotations_df)
        total_all = len(self.genes_df)

        for category in self.functional_categories.keys():
            in_sig = sig_counts.get(category, 0)
            out_sig = total_sig - in_sig
            in_all = all_counts.get(category, 0)
            out_all = total_all - in_all

            # Fisher's exact test
            # [[in_sig, out_sig], [in_all, out_all]]
            odds_ratio, p_value = fisher_exact([[in_sig, out_sig],
                                               [in_all - in_sig, out_all - out_sig]])

            fold_enrichment = (in_sig / total_sig) / (in_all / total_all) if in_all > 0 else 0

            enrichment_results.append({
                'category': category,
                'n_significant': in_sig,
                'n_total': in_all,
                'fold_enrichment': fold_enrichment,
                'p_value': p_value,
                'odds_ratio': odds_ratio
            })

        # Save enrichment results
        enrichment_df = pd.DataFrame(enrichment_results)
        enrichment_df = enrichment_df.sort_values('p_value')
        enrichment_df.to_csv(self.output_dir / 'functional_enrichment.txt',
                            sep='\t', index=False)

        # Log significant enrichments
        significant = enrichment_df[enrichment_df['p_value'] < 0.05]
        if len(significant) > 0:
            self.logger.info("Significant functional enrichments (p < 0.05):")
            for _, row in significant.iterrows():
                self.logger.info(f"  {row['category']}: "
                               f"{row['fold_enrichment']:.2f}x enrichment, "
                               f"p={row['p_value']:.4f}")
        else:
            self.logger.info("No significant functional enrichments found")

    def prioritize_genes(self):
        """Create prioritized candidate gene list."""
        self.logger.info("Prioritizing candidate genes...")

        annotations_df = pd.DataFrame(self.snp_annotations)

        # Define priority levels
        def assign_priority(row):
            """Assign priority based on multiple criteria."""
            p_val = row.get('p', 1.0)
            effect = row.get('effect_size', 1.0)
            category = row['functional_category']
            distance = row['distance']

            # Priority 1: Known stress genes with strong signal
            if category in ['Heat Shock & Stress Response', 'Photosynthesis',
                           'Antioxidant Defense'] and p_val < 1e-8:
                return 1

            # Priority 2: Large effect size
            if abs(effect) > 3.0 or (isinstance(effect, float) and effect > 3.0):
                return 2

            # Priority 3: Close to gene with moderate signal
            if distance < 10000 and p_val < 1e-6:
                return 3

            # Priority 4: Other significant hits
            return 4

        annotations_df['priority'] = annotations_df.apply(assign_priority, axis=1)

        # Sort by priority
        prioritized = annotations_df.sort_values(['priority', 'p'])

        # Save prioritized list
        prioritized.to_csv(self.output_dir / 'priority_genes.csv', index=False)

        # Summary by priority
        priority_summary = prioritized['priority'].value_counts().sort_index()

        summary_file = self.output_dir / 'candidate_genes_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("CANDIDATE GENES SUMMARY\n")
            f.write("="*80 + "\n\n")

            f.write("Total significant SNPs: {}\n\n".format(len(annotations_df)))

            f.write("Priority breakdown:\n")
            for priority, count in priority_summary.items():
                f.write(f"  Priority {priority}: {count} SNPs\n")

            f.write("\nTop 10 candidate genes:\n")
            f.write("-"*80 + "\n")
            for i, (_, row) in enumerate(prioritized.head(10).iterrows(), 1):
                f.write(f"{i}. {row['gene_name']} ({row['gene_id']})\n")
                f.write(f"   Location: {row['chr']}:{row['position']}\n")
                f.write(f"   Category: {row['functional_category']}\n")
                f.write(f"   P-value: {row.get('p', 'N/A')}\n")
                f.write(f"   Description: {row['description']}\n\n")

        self.logger.info(f"Prioritization complete. Top genes saved to {summary_file}")

    def create_visualizations(self):
        """Create annotation visualizations."""
        self.logger.info("Creating visualizations...")

        annotations_df = pd.DataFrame(self.snp_annotations)

        # Set style
        sns.set_style("whitegrid")
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # 1. SNPs by location type
        location_counts = annotations_df['location_type'].value_counts()
        colors = sns.color_palette("Set2", len(location_counts))
        axes[0].pie(location_counts.values, labels=location_counts.index,
                   autopct='%1.1f%%', colors=colors, startangle=90)
        axes[0].set_title('SNPs by Genomic Location')

        # 2. Genes by functional category
        category_counts = annotations_df['functional_category'].value_counts()
        axes[1].barh(range(len(category_counts)), category_counts.values)
        axes[1].set_yticks(range(len(category_counts)))
        axes[1].set_yticklabels(category_counts.index, fontsize=9)
        axes[1].set_xlabel('Number of SNPs')
        axes[1].set_title('SNPs by Functional Category')
        axes[1].invert_yaxis()

        # 3. Enrichment fold-change
        enrichment_df = pd.read_csv(self.output_dir / 'functional_enrichment.txt', sep='\t')
        enrichment_df = enrichment_df[enrichment_df['category'] != 'Unknown/Other']

        # Color by significance
        colors = ['red' if p < 0.05 else 'gray' for p in enrichment_df['p_value']]

        axes[2].barh(range(len(enrichment_df)), enrichment_df['fold_enrichment'].values,
                    color=colors)
        axes[2].set_yticks(range(len(enrichment_df)))
        axes[2].set_yticklabels(enrichment_df['category'], fontsize=9)
        axes[2].set_xlabel('Fold Enrichment')
        axes[2].set_title('Functional Enrichment')
        axes[2].axvline(1.0, color='black', linestyle='--', linewidth=1)
        axes[2].invert_yaxis()

        plt.tight_layout()
        plt.savefig(self.output_dir / 'annotation_plots.png', dpi=300, bbox_inches='tight')
        plt.close()

        self.logger.info(f"Visualizations saved to {self.output_dir / 'annotation_plots.png'}")

    def save_results(self):
        """Save all annotation results."""
        self.logger.info("Saving results...")

        # Save detailed annotations
        annotations_df = pd.DataFrame(self.snp_annotations)
        annotations_df.to_csv(self.output_dir / 'snp_annotations.csv', index=False)

        # Save genes by category
        category_file = self.output_dir / 'genes_by_category.txt'
        with open(category_file, 'w') as f:
            for category in self.functional_categories.keys():
                genes_in_cat = annotations_df[
                    annotations_df['functional_category'] == category
                ]
                f.write(f"\n{category} ({len(genes_in_cat)} SNPs):\n")
                f.write("-"*80 + "\n")
                for _, row in genes_in_cat.iterrows():
                    f.write(f"  {row['gene_name']}: {row['description']}\n")

        self.logger.info(f"Results saved to {self.output_dir}")

    def run_pipeline(self):
        """Run the complete annotation pipeline."""
        self.logger.info("="*80)
        self.logger.info("Starting Gene Annotation Pipeline")
        self.logger.info(f"SNPs file: {self.snps_file}")
        self.logger.info(f"Annotation: {self.annotation_file}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info("="*80)

        # Step 1: Load data
        self.load_snps()
        self.load_annotation()

        # Step 2: Annotate SNPs
        self.annotate_snps()

        # Step 3: Test enrichment
        self.test_enrichment()

        # Step 4: Prioritize genes
        self.prioritize_genes()

        # Step 5: Create visualizations
        self.create_visualizations()

        # Step 6: Save results
        self.save_results()

        self.logger.info("="*80)
        self.logger.info("Gene Annotation Pipeline Completed!")
        self.logger.info(f"Annotated {len(self.snp_annotations)} SNPs")
        self.logger.info(f"Results saved to: {self.output_dir}")
        self.logger.info("="*80)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Gene Annotation Pipeline for Kelp GWAS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python annotate_genes.py --snps gwas_significant.txt \\
                            --annotation Macpyr2_annotation.gff \\
                            --output-dir annotation/ \\
                            --window 50000
        """
    )

    parser.add_argument('--snps', required=True,
                       help='File with significant SNPs (TSV/CSV)')
    parser.add_argument('--annotation', required=True,
                       help='GFF annotation file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('--window', type=int, default=50000,
                       help='Window size for nearby genes (bp, default: 50000)')

    args = parser.parse_args()

    # Create and run pipeline
    pipeline = GeneAnnotator(
        snps_file=args.snps,
        annotation_file=args.annotation,
        output_dir=args.output_dir,
        window_size=args.window
    )

    pipeline.run_pipeline()


if __name__ == '__main__':
    main()
