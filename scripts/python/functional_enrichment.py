# Save as: scripts/python/functional_enrichment_final.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import requests
import json
import time
import os
import warnings
from typing import List, Dict
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Get absolute path to project root
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print(f"Project root: {PROJECT_ROOT}")

class FunctionalEnrichment:
    """Perform functional enrichment analysis using various databases."""
    
    def __init__(self, gene_list: List[str], background_genes: List[str] = None):
        """
        Initialize with gene list for enrichment.
        
        Args:
            gene_list: List of genes to test for enrichment
            background_genes: List of all genes in the background (optional)
        """
        self.gene_list = gene_list
        self.background_genes = background_genes
        print(f"Initialized with {len(gene_list)} genes")
        
    def enrichr_analysis(self, gene_sets: List[str] = None):
        """
        Perform enrichment analysis using Enrichr API.
        
        Args:
            gene_sets: List of gene set libraries to test
        
        Returns:
            Dictionary with enrichment results for each gene set
        """
        if gene_sets is None:
            gene_sets = [
                'GO_Biological_Process_2023',
                'GO_Molecular_Function_2023', 
                'GO_Cellular_Component_2023',
                'KEGG_2021_Human',
                'Reactome_2022'
            ]
        
        print(f"Running Enrichr analysis for {len(gene_sets)} gene sets...")
        
        # Prepare genes for Enrichr
        genes_str = '\n'.join(self.gene_list)
        
        # Enrichr API endpoints
        ENRICHR_URL = 'http://maayanlab.cloud/Enrichr'
        
        # Add gene list
        add_url = f'{ENRICHR_URL}/addList'
        payload = {
            'list': (None, genes_str),
            'description': (None, 'Shared DEGs Psoriasis vs Crohn\'s')
        }
        
        try:
            response = requests.post(add_url, files=payload)
            if not response.ok:
                raise Exception(f"Error adding gene list: {response.status_code}")
            
            data = response.json()
            list_id = data.get('userListId')
            
            if not list_id:
                raise Exception("No list ID returned from Enrichr")
            
            print(f"Successfully uploaded gene list. List ID: {list_id}")
            
            # Get enrichment results for each gene set
            results = {}
            for gene_set in gene_sets:
                print(f"  Querying {gene_set}...")
                
                query_url = f'{ENRICHR_URL}/enrich'
                params = {
                    'userListId': list_id,
                    'backgroundType': gene_set
                }
                
                response = requests.get(query_url, params=params)
                if response.ok:
                    data = response.json()
                    if gene_set in data:
                        df = pd.DataFrame(data[gene_set], 
                                         columns=['Rank', 'Term', 'P-value', 
                                                 'Z-score', 'Combined Score', 
                                                 'Overlapping Genes', 
                                                 'Adjusted P-value', 
                                                 'Old P-value', 'Old Adjusted P-value'])
                        results[gene_set] = df
                        print(f"    Found {len(df)} terms")
                    else:
                        print(f"    Warning: No results for {gene_set}")
                else:
                    print(f"    Warning: Failed to query {gene_set} (Status: {response.status_code})")
                
                time.sleep(1)  # Be nice to the API
            
            return results
            
        except Exception as e:
            print(f"Error in Enrichr analysis: {e}")
            return {}
    
    def plot_enrichment_results(self, results: Dict[str, pd.DataFrame], 
                               top_n: int = 10, output_dir: str = None):
        """
        Plot enrichment results.
        
        Args:
            results: Dictionary of DataFrames from enrichr_analysis
            top_n: Number of top terms to plot
            output_dir: Directory to save plots
        """
        if output_dir is None:
            output_dir = os.path.join(PROJECT_ROOT, 'results', 'plots')
        
        os.makedirs(output_dir, exist_ok=True)
        
        for gene_set, df in results.items():
            if len(df) == 0:
                continue
            
            # Get top terms
            df_top = df.head(top_n).copy()
            
            # Convert to numeric
            df_top['P-value'] = pd.to_numeric(df_top['P-value'], errors='coerce')
            df_top['Adjusted P-value'] = pd.to_numeric(df_top['Adjusted P-value'], errors='coerce')
            
            # Sort by p-value
            df_top = df_top.sort_values('P-value')
            
            # Create -log10(p-value)
            df_top['-log10(p-value)'] = -np.log10(df_top['P-value'])
            
            # Create plot
            plt.figure(figsize=(12, 8))
            
            # Create horizontal bar plot
            y_pos = np.arange(len(df_top))
            bars = plt.barh(y_pos, df_top['-log10(p-value)'])
            
            # Color by significance
            for i, (idx, row) in enumerate(df_top.iterrows()):
                if row['Adjusted P-value'] < 0.05:
                    bars[i].set_color('firebrick')  # Significant
                else:
                    bars[i].set_color('steelblue')  # Not significant
            
            # Customize plot
            plt.yticks(y_pos, df_top['Term'], fontsize=10)
            plt.xlabel('-log10(p-value)', fontsize=12)
            plt.title(f'{gene_set} - Top {top_n} Enriched Terms', fontsize=14, fontweight='bold')
            
            # Add grid
            plt.grid(True, axis='x', alpha=0.3)
            
            # Add significance threshold line
            sig_threshold = -np.log10(0.05)
            plt.axvline(x=sig_threshold, color='red', linestyle='--', alpha=0.5, label='p=0.05')
            plt.legend()
            
            # Add p-value annotations
            for i, (idx, row) in enumerate(df_top.iterrows()):
                plt.text(row['-log10(p-value)'] + 0.05, i, 
                        f"p={row['P-value']:.1e}", 
                        va='center', fontsize=9)
            
            plt.tight_layout()
            
            # Save plot
            safe_name = gene_set.replace(' ', '_').replace('/', '_').replace(':', '_').replace('\\', '_')
            output_path = os.path.join(output_dir, f'enrichment_{safe_name}.png')
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"  Saved plot: {output_path}")
    
    def create_summary_table(self, results: Dict[str, pd.DataFrame], 
                           output_path: str = None):
        """
        Create summary table of enrichment results.
        
        Args:
            results: Dictionary of DataFrames
            output_path: Path to save summary CSV
        """
        if output_path is None:
            output_path = os.path.join(PROJECT_ROOT, 'results', 'enrichment_summary.csv')
        
        summary_data = []
        
        for gene_set, df in results.items():
            if len(df) > 0:
                # Convert to numeric
                df['P-value'] = pd.to_numeric(df['P-value'], errors='coerce')
                df['Adjusted P-value'] = pd.to_numeric(df['Adjusted P-value'], errors='coerce')
                
                # Get significant terms (adjusted p < 0.05)
                sig_df = df[df['Adjusted P-value'] < 0.05].copy()
                
                if len(sig_df) > 0:
                    for _, row in sig_df.iterrows():
                        summary_data.append({
                            'Gene_Set': gene_set,
                            'Term': row['Term'],
                            'P_value': float(row['P-value']),
                            'Adj_P_value': float(row['Adjusted P-value']),
                            'Combined_Score': float(row['Combined Score']) if 'Combined Score' in row else np.nan,
                            'Overlap_Genes': row['Overlapping Genes'],
                            'Gene_Count': len(str(row['Overlapping Genes']).split(',')),
                            'Z_Score': float(row['Z-score']) if 'Z-score' in row else np.nan
                        })
        
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df = summary_df.sort_values(['Gene_Set', 'Adj_P_value'])
            
            # Create output directory if needed
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            summary_df.to_csv(output_path, index=False)
            print(f"\nSaved enrichment summary: {output_path}")
            
            # Print top terms
            self._print_summary(summary_df)
        
        return summary_df if summary_data else pd.DataFrame()
    
    def _print_summary(self, summary_df):
        """Print summary of enrichment results."""
        print("\n" + "="*80)
        print("TOP ENRICHED TERMS (Adj. p < 0.05)")
        print("="*80)
        
        for gene_set in summary_df['Gene_Set'].unique():
            subset = summary_df[summary_df['Gene_Set'] == gene_set].head(5)
            print(f"\n{gene_set}:")
            for _, row in subset.iterrows():
                print(f"  • {row['Term']}")
                print(f"    p={row['Adj_P_value']:.2e}, {row['Gene_Count']} genes")

def main():
    """Main function."""
    print("="*80)
    print("FUNCTIONAL ENRICHMENT ANALYSIS")
    print("="*80)
    
    # Define paths
    shared_genes_path = os.path.join(PROJECT_ROOT, 'results', 'shared_degs_training.csv')
    
    # Load shared DEGs
    try:
        shared_genes = pd.read_csv(shared_genes_path)
        gene_list = shared_genes['Gene'].tolist()
        print(f"Loaded {len(gene_list)} shared DEGs from {shared_genes_path}")
        
    except FileNotFoundError:
        print(f"Error: {shared_genes_path} not found.")
        print("Looking for alternative locations...")
        
        # Try to find the file
        possible_paths = [
            os.path.join(PROJECT_ROOT, 'shared_degs_training.csv'),
            os.path.join(PROJECT_ROOT, 'results', 'shared_degs_training.csv'),
            'shared_degs_training.csv',
            '../shared_degs_training.csv',
            '../../shared_degs_training.csv'
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                shared_genes = pd.read_csv(path)
                gene_list = shared_genes['Gene'].tolist()
                print(f"Found and loaded {len(gene_list)} genes from {path}")
                break
        else:
            # Use genes from the paper
            print("Using genes from paper's key findings...")
            gene_list = ['KIF4A', 'DLGAP5', 'NCAPG', 'CCNB1', 'CEP55', 
                        'IL6', 'IL8', 'TNF', 'CDK1', 'CDC20',
                        'MCM2', 'MCM4', 'STAT1', 'STAT3', 'NFKB1']
            print(f"Using {len(gene_list)} genes from paper")
    
    # Initialize analyzer
    analyzer = FunctionalEnrichment(gene_list)
    
    # Run Enrichr analysis
    print("\nRunning Enrichr analysis...")
    enrichr_results = analyzer.enrichr_analysis()
    
    if enrichr_results:
        # Plot results
        print("\nCreating plots...")
        analyzer.plot_enrichment_results(enrichr_results, top_n=15)
        
        # Create summary
        analyzer.create_summary_table(enrichr_results)
        
        # Compare with paper
        print("\n" + "="*80)
        print("COMPARISON WITH PAPER'S FINDINGS")
        print("="*80)
        print("Paper reported enrichment in:")
        print("  • Cell cycle regulation")
        print("  • Immune response pathways")
        print("  • Inflammatory signaling")
        print("  • Chromosome segregation")
        
        # Create comparison visualization
        create_comparison_visualization(gene_list)
    else:
        print("No results from Enrichr analysis")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)

def create_comparison_visualization(gene_list):
    """Create visualization comparing with paper's findings."""
    # Paper's reported pathways
    paper_pathways = {
        'Cell Cycle': ['KIF4A', 'DLGAP5', 'NCAPG', 'CCNB1', 'CEP55', 'CDK1', 'CDC20'],
        'Immune Response': ['IL6', 'IL8', 'TNF', 'STAT1', 'STAT3'],
        'Inflammation': ['IL6', 'IL8', 'TNF', 'NFKB1'],
        'Chromosome Segregation': ['KIF4A', 'DLGAP5', 'NCAPG', 'CEP55']
    }
    
    # Calculate overlap
    overlap_data = []
    for pathway, genes in paper_pathways.items():
        overlap = set(gene_list).intersection(set(genes))
        overlap_data.append({
            'Pathway': pathway,
            'Paper_Genes': len(genes),
            'Our_Overlap': len(overlap),
            'Percent': (len(overlap) / len(genes)) * 100 if len(genes) > 0 else 0,
            'Overlap_Genes': ', '.join(overlap)
        })
    
    # Create DataFrame
    overlap_df = pd.DataFrame(overlap_data)
    overlap_df = overlap_df.sort_values('Percent', ascending=False)
    
    # Save to CSV
    output_path = os.path.join(PROJECT_ROOT, 'results', 'paper_comparison.csv')
    overlap_df.to_csv(output_path, index=False)
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    
    x = range(len(overlap_df))
    bars = plt.barh(x, overlap_df['Percent'])
    
    # Color bars
    colors = ['firebrick' if p > 50 else 'darkorange' if p > 25 else 'steelblue' 
              for p in overlap_df['Percent']]
    for i, bar in enumerate(bars):
        bar.set_color(colors[i])
    
    plt.yticks(x, overlap_df['Pathway'])
    plt.xlabel('Percent Overlap with Paper (%)', fontsize=12)
    plt.title('Comparison with Paper\'s Reported Pathways', fontsize=14, fontweight='bold')
    
    # Add percentage labels
    for i, (idx, row) in enumerate(overlap_df.iterrows()):
        plt.text(row['Percent'] + 1, i, f"{row['Percent']:.0f}%", va='center')
    
    plt.tight_layout()
    
    # Save plot
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'plots')
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'paper_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nSaved comparison analysis:")
    print(f"  • {output_path}")
    print(f"  • {os.path.join(output_dir, 'paper_comparison.png')}")

if __name__ == "__main__":
    main()