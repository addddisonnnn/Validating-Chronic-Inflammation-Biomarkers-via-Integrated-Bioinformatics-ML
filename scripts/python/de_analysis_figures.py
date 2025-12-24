# Save as: scripts/python/de_analysis_figures.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
from statsmodels.stats.multitest import multipletests
from matplotlib_venn import venn2
import os
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Get project root directory
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print(f"Project root: {PROJECT_ROOT}")

def setup_output_directory():
    """Create output directory for all figures."""
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'Figure2_DE_Analysis')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")
    return output_dir

def find_file(filename, search_dirs=None):
    """Find file in various possible locations."""
    if search_dirs is None:
        search_dirs = [
            os.path.join(PROJECT_ROOT, 'data', 'processed'),
            os.path.join(PROJECT_ROOT, 'data'),
            os.path.join(PROJECT_ROOT, 'results'),
            os.path.join(PROJECT_ROOT, 'metadata'),
            PROJECT_ROOT
        ]
    
    for search_dir in search_dirs:
        path = os.path.join(search_dir, filename)
        if os.path.exists(path):
            return path
    
    # Also check recursively
    for root, dirs, files in os.walk(PROJECT_ROOT):
        if filename in files:
            return os.path.join(root, filename)
    
    return None

def load_expression_data(dataset_name):
    """Load expression data for a dataset."""
    print(f"\nLoading {dataset_name} expression data...")
    
    # Try different filename patterns
    patterns = [
        f"{dataset_name}_normalized_gene_expression.csv",
        f"{dataset_name}_expression.csv",
        f"{dataset_name}.csv"
    ]
    
    for pattern in patterns:
        filepath = find_file(pattern)
        if filepath:
            try:
                expr_data = pd.read_csv(filepath, index_col=0)
                print(f"  ✓ Loaded from: {filepath}")
                print(f"    Shape: {expr_data.shape[0]} genes × {expr_data.shape[1]} samples")
                return expr_data
            except Exception as e:
                print(f"  ✗ Error loading {filepath}: {e}")
    
    print(f"  ✗ Could not find expression data for {dataset_name}")
    return None

def load_deg_results(dataset_name):
    """Load DEG results for a dataset."""
    print(f"Loading {dataset_name} DEG results...")
    
    patterns = [
        f"{dataset_name}_DEGs.csv",
        f"{dataset_name}_DEGs_python.csv",
        f"{dataset_name}_degs.csv"
    ]
    
    for pattern in patterns:
        filepath = find_file(pattern)
        if filepath:
            try:
                degs = pd.read_csv(filepath)
                print(f"  ✓ Loaded from: {filepath}")
                print(f"    Found {len(degs)} DEGs")
                return degs
            except Exception as e:
                print(f"  ✗ Error loading {filepath}: {e}")
    
    print(f"  ✗ Could not find DEG results for {dataset_name}")
    return None

def load_sample_metadata(dataset_name):
    """Load sample metadata."""
    print(f"Loading {dataset_name} metadata...")
    
    if 'GSE13355' in dataset_name:
        patterns = ['GSE13355_sample_groups.csv']
        disease_label = 'Psoriasis'
    elif 'GSE75214' in dataset_name:
        patterns = ['GSE75214_sample_groups_fixed.csv', 'GSE75214_sample_groups.csv']
        disease_label = 'CD'
    else:
        print(f"  ✗ Unknown dataset: {dataset_name}")
        return None, [], []
    
    for pattern in patterns:
        filepath = find_file(pattern)
        if filepath:
            try:
                metadata = pd.read_csv(filepath)
                print(f"  ✓ Loaded from: {filepath}")
                
                # Get disease and control samples
                disease_samples = metadata[metadata['Group'] == disease_label]['GSM_ID'].tolist()
                control_samples = metadata[metadata['Group'] == 'Control']['GSM_ID'].tolist()
                
                print(f"    Disease: {len(disease_samples)} samples")
                print(f"    Control: {len(control_samples)} samples")
                
                return metadata, disease_samples, control_samples
            except Exception as e:
                print(f"  ✗ Error loading {filepath}: {e}")
    
    print(f"  ✗ Could not find metadata for {dataset_name}")
    return None, [], []

def match_samples_to_expression(expr_data, disease_samples, control_samples):
    """Match sample IDs from metadata to expression data columns."""
    expr_columns = expr_data.columns.tolist()
    
    # Function to find matching column
    def find_matching_column(sample_id):
        # Exact match
        if sample_id in expr_columns:
            return sample_id
        
        # Partial match (sample_id might be part of column name)
        for col in expr_columns:
            if sample_id in col:
                return col
        
        # Try without prefix/suffix
        clean_sample = sample_id.replace('GSM', '').strip()
        for col in expr_columns:
            if clean_sample in col:
                return col
        
        return None
    
    # Match disease samples
    matched_disease = []
    for sample in disease_samples:
        match = find_matching_column(sample)
        if match:
            matched_disease.append(match)
    
    # Match control samples
    matched_control = []
    for sample in control_samples:
        match = find_matching_column(sample)
        if match:
            matched_control.append(match)
    
    print(f"  Matched {len(matched_disease)}/{len(disease_samples)} disease samples")
    print(f"  Matched {len(matched_control)}/{len(control_samples)} control samples")
    
    return matched_disease, matched_control

def create_pca_figure(psoriasis_expr, psoriasis_disease, psoriasis_control,
                     crohns_expr, crohns_disease, crohns_control, output_dir):
    """Create PCA plots (Figure 2A and 2B)."""
    print("\nCreating PCA plots...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot A: Psoriasis PCA
    if len(psoriasis_disease) > 1 and len(psoriasis_control) > 1:
        try:
            # Prepare data
            samples = psoriasis_disease + psoriasis_control
            X = psoriasis_expr[samples].T.values
            X_scaled = StandardScaler().fit_transform(X)
            
            # PCA
            pca = PCA(n_components=2)
            pc = pca.fit_transform(X_scaled)
            
            # Create plot
            scatter1 = axes[0].scatter(pc[:len(psoriasis_disease), 0], 
                                      pc[:len(psoriasis_disease), 1],
                                      c='red', label='Psoriasis', alpha=0.7, s=50)
            scatter2 = axes[0].scatter(pc[len(psoriasis_disease):, 0],
                                      pc[len(psoriasis_disease):, 1],
                                      c='blue', label='Control', alpha=0.7, s=50)
            
            axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
            axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
            axes[0].legend()
            axes[0].set_title('A) PCA: Psoriasis (GSE13355)', fontweight='bold', fontsize=12)
            axes[0].grid(True, alpha=0.3)
            
        except Exception as e:
            axes[0].text(0.5, 0.5, f'Error: {str(e)[:50]}', 
                        ha='center', va='center', transform=axes[0].transAxes)
            axes[0].set_title('A) PCA: Psoriasis (Error)', fontweight='bold', fontsize=12)
    else:
        axes[0].text(0.5, 0.5, 'Insufficient data', 
                    ha='center', va='center', transform=axes[0].transAxes)
        axes[0].set_title('A) PCA: Psoriasis', fontweight='bold', fontsize=12)
    
    # Plot B: Crohn's PCA
    if len(crohns_disease) > 1 and len(crohns_control) > 1:
        try:
            # Prepare data
            samples = crohns_disease + crohns_control
            X = crohns_expr[samples].T.values
            X_scaled = StandardScaler().fit_transform(X)
            
            # PCA
            pca = PCA(n_components=2)
            pc = pca.fit_transform(X_scaled)
            
            # Create plot
            axes[1].scatter(pc[:len(crohns_disease), 0], 
                          pc[:len(crohns_disease), 1],
                          c='red', label='Crohn\'s Disease', alpha=0.7, s=50)
            axes[1].scatter(pc[len(crohns_disease):, 0],
                          pc[len(crohns_disease):, 1],
                          c='blue', label='Control', alpha=0.7, s=50)
            
            axes[1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
            axes[1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
            axes[1].legend()
            axes[1].set_title('B) PCA: Crohn\'s Disease (GSE75214)', fontweight='bold', fontsize=12)
            axes[1].grid(True, alpha=0.3)
            
        except Exception as e:
            axes[1].text(0.5, 0.5, f'Error: {str(e)[:50]}', 
                        ha='center', va='center', transform=axes[1].transAxes)
            axes[1].set_title('B) PCA: Crohn\'s (Error)', fontweight='bold', fontsize=12)
    else:
        axes[1].text(0.5, 0.5, 'Insufficient data', 
                    ha='center', va='center', transform=axes[1].transAxes)
        axes[1].set_title('B) PCA: Crohn\'s Disease', fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure2_AB_PCA.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_path}")
    plt.close()
    
    # Also save individual plots
    for i, (title, suffix) in enumerate([('Psoriasis', 'A'), ('Crohns', 'B')]):
        fig, ax = plt.subplots(figsize=(8, 6))
        if i == 0 and len(psoriasis_disease) > 1:
            ax.scatter(pc[:len(psoriasis_disease), 0], pc[:len(psoriasis_disease), 1],
                      c='red', label='Disease', alpha=0.7, s=50)
            ax.scatter(pc[len(psoriasis_disease):, 0], pc[len(psoriasis_disease):, 1],
                      c='blue', label='Control', alpha=0.7, s=50)
            ax.set_title(f'PCA: {title}', fontweight='bold')
        elif i == 1 and len(crohns_disease) > 1:
            ax.scatter(pc[:len(crohns_disease), 0], pc[:len(crohns_disease), 1],
                      c='red', label='Disease', alpha=0.7, s=50)
            ax.scatter(pc[len(crohns_disease):, 0], pc[len(crohns_disease):, 1],
                      c='blue', label='Control', alpha=0.7, s=50)
            ax.set_title(f'PCA: {title}', fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'PCA: {title}', fontweight='bold')
        
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        individual_path = os.path.join(output_dir, f'Figure2{suffix}_PCA_{title}.png')
        plt.savefig(individual_path, dpi=300, bbox_inches='tight')
        plt.close()

def create_volcano_figure(psoriasis_degs, crohns_degs, output_dir):
    """Create volcano plots (Figure 2C and 2D)."""
    print("\nCreating volcano plots...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot C: Psoriasis Volcano
    if psoriasis_degs is not None and len(psoriasis_degs) > 0:
        try:
            # Check if we have necessary columns
            if 'log2FC' in psoriasis_degs.columns and 'Adj_PValue' in psoriasis_degs.columns:
                data = psoriasis_degs.copy()
                data['-log10(p)'] = -np.log10(data['Adj_PValue'])
                data['Significant'] = (abs(data['log2FC']) > 0.585) & (data['Adj_PValue'] < 0.05)
                data['Color'] = np.where(data['Significant'], 
                                        np.where(data['log2FC'] > 0, 'red', 'blue'), 
                                        'gray')
                
                # Plot
                scatter = axes[0].scatter(data['log2FC'], data['-log10(p)'],
                                         c=data['Color'], alpha=0.6, s=30)
                
                # Add lines
                axes[0].axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, linewidth=1)
                axes[0].axvline(x=0.585, color='black', linestyle='--', alpha=0.5, linewidth=1)
                axes[0].axvline(x=-0.585, color='black', linestyle='--', alpha=0.5, linewidth=1)
                
                # Add legend manually
                from matplotlib.patches import Patch
                legend_elements = [
                    Patch(facecolor='red', alpha=0.6, label='Upregulated'),
                    Patch(facecolor='blue', alpha=0.6, label='Downregulated'),
                    Patch(facecolor='gray', alpha=0.6, label='Not significant')
                ]
                axes[0].legend(handles=legend_elements, loc='upper right')
                
                axes[0].set_xlabel('log2 Fold Change')
                axes[0].set_ylabel('-log10(Adjusted P-value)')
                axes[0].set_title('C) Volcano: Psoriasis DEGs', fontweight='bold', fontsize=12)
                axes[0].grid(True, alpha=0.3)
                
                # Count significant
                sig_up = sum((data['log2FC'] > 0.585) & (data['Adj_PValue'] < 0.05))
                sig_down = sum((data['log2FC'] < -0.585) & (data['Adj_PValue'] < 0.05))
                axes[0].text(0.02, 0.98, f'Up: {sig_up}\nDown: {sig_down}',
                            transform=axes[0].transAxes, fontsize=10,
                            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            else:
                axes[0].text(0.5, 0.5, 'Missing log2FC or P-value columns',
                            ha='center', va='center', transform=axes[0].transAxes)
                axes[0].set_title('C) Volcano: Psoriasis', fontweight='bold', fontsize=12)
                
        except Exception as e:
            axes[0].text(0.5, 0.5, f'Error: {str(e)[:50]}',
                        ha='center', va='center', transform=axes[0].transAxes)
            axes[0].set_title('C) Volcano: Psoriasis (Error)', fontweight='bold', fontsize=12)
    else:
        axes[0].text(0.5, 0.5, 'No DEG data available',
                    ha='center', va='center', transform=axes[0].transAxes)
        axes[0].set_title('C) Volcano: Psoriasis', fontweight='bold', fontsize=12)
    
    # Plot D: Crohn's Volcano
    if crohns_degs is not None and len(crohns_degs) > 0:
        try:
            if 'log2FC' in crohns_degs.columns and 'Adj_PValue' in crohns_degs.columns:
                data = crohns_degs.copy()
                data['-log10(p)'] = -np.log10(data['Adj_PValue'])
                data['Significant'] = (abs(data['log2FC']) > 0.585) & (data['Adj_PValue'] < 0.05)
                data['Color'] = np.where(data['Significant'],
                                        np.where(data['log2FC'] > 0, 'red', 'blue'),
                                        'gray')
                
                # Plot
                axes[1].scatter(data['log2FC'], data['-log10(p)'],
                               c=data['Color'], alpha=0.6, s=30)
                
                # Add lines
                axes[1].axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, linewidth=1)
                axes[1].axvline(x=0.585, color='black', linestyle='--', alpha=0.5, linewidth=1)
                axes[1].axvline(x=-0.585, color='black', linestyle='--', alpha=0.5, linewidth=1)
                
                # Legend
                from matplotlib.patches import Patch
                legend_elements = [
                    Patch(facecolor='red', alpha=0.6, label='Upregulated'),
                    Patch(facecolor='blue', alpha=0.6, label='Downregulated'),
                    Patch(facecolor='gray', alpha=0.6, label='Not significant')
                ]
                axes[1].legend(handles=legend_elements, loc='upper right')
                
                axes[1].set_xlabel('log2 Fold Change')
                axes[1].set_ylabel('-log10(Adjusted P-value)')
                axes[1].set_title('D) Volcano: Crohn\'s Disease DEGs', fontweight='bold', fontsize=12)
                axes[1].grid(True, alpha=0.3)
                
                # Count significant
                sig_up = sum((data['log2FC'] > 0.585) & (data['Adj_PValue'] < 0.05))
                sig_down = sum((data['log2FC'] < -0.585) & (data['Adj_PValue'] < 0.05))
                axes[1].text(0.02, 0.98, f'Up: {sig_up}\nDown: {sig_down}',
                            transform=axes[1].transAxes, fontsize=10,
                            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            else:
                axes[1].text(0.5, 0.5, 'Missing log2FC or P-value columns',
                            ha='center', va='center', transform=axes[1].transAxes)
                axes[1].set_title('D) Volcano: Crohn\'s', fontweight='bold', fontsize=12)
                
        except Exception as e:
            axes[1].text(0.5, 0.5, f'Error: {str(e)[:50]}',
                        ha='center', va='center', transform=axes[1].transAxes)
            axes[1].set_title('D) Volcano: Crohn\'s (Error)', fontweight='bold', fontsize=12)
    else:
        axes[1].text(0.5, 0.5, 'No DEG data available',
                    ha='center', va='center', transform=axes[1].transAxes)
        axes[1].set_title('D) Volcano: Crohn\'s Disease', fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure2_CD_Volcano.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_path}")
    plt.close()

def create_venn_diagram_figure(psoriasis_degs, crohns_degs, output_dir):
    """Create Venn diagram (Figure 2E)."""
    print("\nCreating Venn diagram...")
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    if psoriasis_degs is not None and crohns_degs is not None:
        try:
            # Get gene sets
            psoriasis_genes = set(psoriasis_degs['Gene'].astype(str).str.upper())
            crohns_genes = set(crohns_degs['Gene'].astype(str).str.upper())
            
            # Create Venn diagram
            v = venn2([psoriasis_genes, crohns_genes],
                     set_labels=('Psoriasis', 'Crohn\'s Disease'), ax=ax)
            
            # Customize colors
            if v.get_patch_by_id('10'):
                v.get_patch_by_id('10').set_color('lightcoral')
                v.get_patch_by_id('10').set_alpha(0.7)
            if v.get_patch_by_id('01'):
                v.get_patch_by_id('01').set_color('lightblue')
                v.get_patch_by_id('01').set_alpha(0.7)
            if v.get_patch_by_id('11'):
                v.get_patch_by_id('11').set_color('lightgreen')
                v.get_patch_by_id('11').set_alpha(0.7)
            
            # Update labels with counts
            v.get_label_by_id('10').set_text(f'Psoriasis only\n{len(psoriasis_genes - crohns_genes)}')
            v.get_label_by_id('01').set_text(f'CD only\n{len(crohns_genes - psoriasis_genes)}')
            v.get_label_by_id('11').set_text(f'Shared\n{len(psoriasis_genes & crohns_genes)}')
            
            # Make labels larger
            for label in v.set_labels:
                if label:
                    label.set_fontsize(12)
                    label.set_fontweight('bold')
            
            for label in v.subset_labels:
                if label:
                    label.set_fontsize(11)
            
            ax.set_title('E) Shared DEGs: Psoriasis vs Crohn\'s Disease', 
                        fontweight='bold', fontsize=14, pad=20)
            
            # Calculate and display shared genes count
            shared_genes = psoriasis_genes & crohns_genes
            print(f"  Shared DEGs: {len(shared_genes)}")
            print(f"  Psoriasis unique: {len(psoriasis_genes - crohns_genes)}")
            print(f"  Crohn's unique: {len(crohns_genes - psoriasis_genes)}")
            
            # Save shared genes list
            if len(shared_genes) > 0:
                shared_df = pd.DataFrame({'Gene': sorted(list(shared_genes))})
                shared_path = os.path.join(output_dir, 'shared_deg_genes.csv')
                shared_df.to_csv(shared_path, index=False)
                print(f"  ✓ Saved shared genes list: {shared_path}")
                
        except Exception as e:
            ax.text(0.5, 0.5, f'Error creating Venn: {str(e)[:50]}',
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title('E) Shared DEGs (Error)', fontweight='bold', fontsize=14)
    else:
        ax.text(0.5, 0.5, 'Missing DEG data for one or both datasets',
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title('E) Shared DEGs', fontweight='bold', fontsize=14)
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure2E_Venn.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_path}")
    plt.close()

def create_heatmap_figures(psoriasis_expr, crohns_expr, 
                          psoriasis_disease, psoriasis_control,
                          crohns_disease, crohns_control,
                          shared_genes, output_dir):
    """Create heatmaps (Figure 2F and 2G)."""
    print("\nCreating heatmaps...")
    
    # Load shared genes if not provided
    if shared_genes is None:
        shared_path = os.path.join(output_dir, 'shared_deg_genes.csv')
        if os.path.exists(shared_path):
            shared_df = pd.read_csv(shared_path)
            shared_genes = shared_df['Gene'].tolist()
        else:
            print("  No shared genes file found")
            shared_genes = []
    
    if len(shared_genes) == 0:
        print("  No shared genes for heatmaps")
        return
    
    # Limit to genes present in both datasets
    psoriasis_genes = [g for g in shared_genes if g in psoriasis_expr.index]
    crohns_genes = [g for g in shared_genes if g in crohns_expr.index]
    
    print(f"  Genes for psoriasis heatmap: {len(psoriasis_genes)}")
    print(f"  Genes for Crohn's heatmap: {len(crohns_genes)}")
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=(14, 12))
    
    # Plot F: Psoriasis heatmap
    if len(psoriasis_genes) > 0 and len(psoriasis_disease) > 0 and len(psoriasis_control) > 0:
        try:
            # Get expression data
            samples = psoriasis_disease + psoriasis_control
            expr_subset = psoriasis_expr.loc[psoriasis_genes, samples]
            
            # Z-score normalize by gene
            expr_z = expr_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
            
            # Create heatmap
            sns.heatmap(expr_z, cmap='RdBu_r', center=0, 
                       yticklabels=False, ax=axes[0],
                       cbar_kws={'label': 'Z-score', 'shrink': 0.8})
            
            # Add separator line
            axes[0].axvline(x=len(psoriasis_disease), color='black', linewidth=2)
            
            # Add labels
            axes[0].set_xlabel('Samples', fontsize=11)
            axes[0].set_ylabel(f'{len(psoriasis_genes)} Shared DEGs', fontsize=11)
            axes[0].set_title('F) Psoriasis: Shared DEGs Expression', 
                            fontweight='bold', fontsize=12, pad=10)
            
            # Add group annotations
            axes[0].text(len(psoriasis_disease)/2, -0.05, 'Psoriasis', 
                        ha='center', va='top', transform=axes[0].get_xaxis_transform(),
                        fontsize=10, fontweight='bold')
            axes[0].text(len(psoriasis_disease) + len(psoriasis_control)/2, -0.05, 'Control',
                        ha='center', va='top', transform=axes[0].get_xaxis_transform(),
                        fontsize=10, fontweight='bold')
            
        except Exception as e:
            axes[0].text(0.5, 0.5, f'Error: {str(e)[:50]}',
                        ha='center', va='center', transform=axes[0].transAxes)
            axes[0].set_title('F) Psoriasis Heatmap (Error)', fontweight='bold', fontsize=12)
    else:
        axes[0].text(0.5, 0.5, 'Insufficient data for heatmap',
                    ha='center', va='center', transform=axes[0].transAxes)
        axes[0].set_title('F) Psoriasis: Shared DEGs', fontweight='bold', fontsize=12)
    
    # Plot G: Crohn's heatmap
    if len(crohns_genes) > 0 and len(crohns_disease) > 0 and len(crohns_control) > 0:
        try:
            # Get expression data
            samples = crohns_disease + crohns_control
            expr_subset = crohns_expr.loc[crohns_genes, samples]
            
            # Z-score normalize by gene
            expr_z = expr_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
            
            # Create heatmap
            sns.heatmap(expr_z, cmap='RdBu_r', center=0,
                       yticklabels=False, ax=axes[1],
                       cbar_kws={'label': 'Z-score', 'shrink': 0.8})
            
            # Add separator line
            axes[1].axvline(x=len(crohns_disease), color='black', linewidth=2)
            
            # Add labels
            axes[1].set_xlabel('Samples', fontsize=11)
            axes[1].set_ylabel(f'{len(crohns_genes)} Shared DEGs', fontsize=11)
            axes[1].set_title('G) Crohn\'s Disease: Shared DEGs Expression',
                            fontweight='bold', fontsize=12, pad=10)
            
            # Add group annotations
            axes[1].text(len(crohns_disease)/2, -0.05, 'Crohn\'s Disease',
                        ha='center', va='top', transform=axes[1].get_xaxis_transform(),
                        fontsize=10, fontweight='bold')
            axes[1].text(len(crohns_disease) + len(crohns_control)/2, -0.05, 'Control',
                        ha='center', va='top', transform=axes[1].get_xaxis_transform(),
                        fontsize=10, fontweight='bold')
            
        except Exception as e:
            axes[1].text(0.5, 0.5, f'Error: {str(e)[:50]}',
                        ha='center', va='center', transform=axes[1].transAxes)
            axes[1].set_title('G) Crohn\'s Heatmap (Error)', fontweight='bold', fontsize=12)
    else:
        axes[1].text(0.5, 0.5, 'Insufficient data for heatmap',
                    ha='center', va='center', transform=axes[1].transAxes)
        axes[1].set_title('G) Crohn\'s: Shared DEGs', fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure2_FG_Heatmaps.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_path}")
    plt.close()

def create_combined_figure(output_dir):
    """Create a combined figure with all plots."""
    print("\nCreating combined figure...")
    
    # This would combine all individual figures into one
    # For now, we'll create a summary figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    # Add titles for each subplot location
    titles = [
        'A) PCA: Psoriasis',
        'B) PCA: Crohn\'s Disease', 
        'C) Volcano: Psoriasis',
        'D) Volcano: Crohn\'s',
        'E) Shared DEGs',
        'F) Expression Heatmaps'
    ]
    
    for i, (ax, title) in enumerate(zip(axes, titles)):
        ax.text(0.5, 0.5, f'{title}\n(See individual plots for details)',
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title(title, fontweight='bold', fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.suptitle('Figure 2: Differential Expression Analysis Results', 
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, 'Figure2_Combined_Summary.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved combined summary: {output_path}")
    plt.close()

def main():
    """Main function to create all figures."""
    print("="*80)
    print("CREATING FIGURE 2: DIFFERENTIAL EXPRESSION ANALYSIS")
    print("="*80)
    
    # Setup output directory
    output_dir = setup_output_directory()
    
    # Load data
    print("\n" + "="*80)
    print("LOADING DATA")
    print("="*80)
    
    # Psoriasis data
    psoriasis_expr = load_expression_data("GSE13355")
    psoriasis_meta, psoriasis_disease_raw, psoriasis_control_raw = load_sample_metadata("GSE13355")
    psoriasis_degs = load_deg_results("GSE13355")
    
    # Crohn's data
    crohns_expr = load_expression_data("GSE75214")
    crohns_meta, crohns_disease_raw, crohns_control_raw = load_sample_metadata("GSE75214")
    crohns_degs = load_deg_results("GSE75214")
    
    # Match samples to expression data
    print("\n" + "="*80)
    print("MATCHING SAMPLES TO EXPRESSION DATA")
    print("="*80)
    
    if psoriasis_expr is not None:
        psoriasis_disease, psoriasis_control = match_samples_to_expression(
            psoriasis_expr, psoriasis_disease_raw, psoriasis_control_raw)
    else:
        psoriasis_disease, psoriasis_control = [], []
    
    if crohns_expr is not None:
        crohns_disease, crohns_control = match_samples_to_expression(
            crohns_expr, crohns_disease_raw, crohns_control_raw)
    else:
        crohns_disease, crohns_control = [], []
    
    # Create figures
    print("\n" + "="*80)
    print("CREATING FIGURES")
    print("="*80)
    
    # Figure 2A-B: PCA plots
    create_pca_figure(psoriasis_expr, psoriasis_disease, psoriasis_control,
                     crohns_expr, crohns_disease, crohns_control, output_dir)
    
    # Figure 2C-D: Volcano plots
    create_volcano_figure(psoriasis_degs, crohns_degs, output_dir)
    
    # Figure 2E: Venn diagram
    create_venn_diagram_figure(psoriasis_degs, crohns_degs, output_dir)
    
    # Figure 2F-G: Heatmaps (need to get shared genes from Venn)
    # We'll create these after getting shared genes
    if psoriasis_degs is not None and crohns_degs is not None:
        psoriasis_genes = set(psoriasis_degs['Gene'].astype(str).str.upper())
        crohns_genes = set(crohns_degs['Gene'].astype(str).str.upper())
        shared_genes = list(psoriasis_genes & crohns_genes)
        
        create_heatmap_figures(psoriasis_expr, crohns_expr,
                              psoriasis_disease, psoriasis_control,
                              crohns_disease, crohns_control,
                              shared_genes, output_dir)
    
    # Create combined summary figure
    create_combined_figure(output_dir)
    
    # Create README file for the output
    create_readme_file(output_dir, 
                      len(psoriasis_disease) if psoriasis_expr is not None else 0,
                      len(psoriasis_control) if psoriasis_expr is not None else 0,
                      len(crohns_disease) if crohns_expr is not None else 0,
                      len(crohns_control) if crohns_expr is not None else 0,
                      psoriasis_degs, crohns_degs)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nAll figures saved to: {output_dir}")
    print("\nFiles created:")
    print("  • Figure2_AB_PCA.png        - PCA plots")
    print("  • Figure2_CD_Volcano.png    - Volcano plots")
    print("  • Figure2E_Venn.png         - Venn diagram")
    print("  • Figure2_FG_Heatmaps.png   - Expression heatmaps")
    print("  • Figure2_Combined_Summary.png - Summary figure")
    print("  • shared_deg_genes.csv      - List of shared DEGs")
    print("  • README.md                 - Analysis summary")

def create_readme_file(output_dir, psoriasis_disease_count, psoriasis_control_count,
                      crohns_disease_count, crohns_control_count,
                      psoriasis_degs, crohns_degs):
    """Create README file with analysis summary."""
    readme_path = os.path.join(output_dir, 'README.md')
    
    with open(readme_path, 'w') as f:
        f.write("# Figure 2: Differential Expression Analysis\n\n")
        f.write("## Summary\n\n")
        
        f.write("### Dataset Information\n")
        f.write(f"- **Psoriasis (GSE13355)**: {psoriasis_disease_count} disease vs {psoriasis_control_count} control samples\n")
        f.write(f"- **Crohn's Disease (GSE75214)**: {crohns_disease_count} disease vs {crohns_control_count} control samples\n\n")
        
        if psoriasis_degs is not None:
            f.write(f"### Psoriasis DEGs: {len(psoriasis_degs)}\n")
            if 'log2FC' in psoriasis_degs.columns:
                up = sum(psoriasis_degs['log2FC'] > 0)
                down = sum(psoriasis_degs['log2FC'] < 0)
                f.write(f"- Upregulated: {up}\n")
                f.write(f"- Downregulated: {down}\n\n")
        
        if crohns_degs is not None:
            f.write(f"### Crohn's Disease DEGs: {len(crohns_degs)}\n")
            if 'log2FC' in crohns_degs.columns:
                up = sum(crohns_degs['log2FC'] > 0)
                down = sum(crohns_degs['log2FC'] < 0)
                f.write(f"- Upregulated: {up}\n")
                f.write(f"- Downregulated: {down}\n\n")
        
        if psoriasis_degs is not None and crohns_degs is not None:
            psoriasis_genes = set(psoriasis_degs['Gene'].astype(str).str.upper())
            crohns_genes = set(crohns_degs['Gene'].astype(str).str.upper())
            shared = len(psoriasis_genes & crohns_genes)
            f.write(f"### Shared DEGs: {shared}\n")
            f.write(f"- Psoriasis unique: {len(psoriasis_genes - crohns_genes)}\n")
            f.write(f"- Crohn's unique: {len(crohns_genes - psoriasis_genes)}\n\n")
        
        f.write("## Files in this Directory\n\n")
        f.write("1. **Figure2_AB_PCA.png** - Principal Component Analysis plots\n")
        f.write("   - A) Psoriasis (GSE13355)\n")
        f.write("   - B) Crohn's Disease (GSE75214)\n\n")
        
        f.write("2. **Figure2_CD_Volcano.png** - Volcano plots of differentially expressed genes\n")
        f.write("   - C) Psoriasis DEGs\n")
        f.write("   - D) Crohn's Disease DEGs\n\n")
        
        f.write("3. **Figure2E_Venn.png** - Venn diagram of shared DEGs\n")
        f.write("   - E) Overlap between psoriasis and Crohn's disease\n\n")
        
        f.write("4. **Figure2_FG_Heatmaps.png** - Expression heatmaps of shared DEGs\n")
        f.write("   - F) Psoriasis samples\n")
        f.write("   - G) Crohn's Disease samples\n\n")
        
        f.write("5. **Figure2_Combined_Summary.png** - Combined overview figure\n\n")
        
        f.write("6. **shared_deg_genes.csv** - List of shared differentially expressed genes\n\n")
        
        f.write("## Analysis Parameters\n\n")
        f.write("- **Differential expression criteria**: |log2FC| > 0.585 and adjusted p-value < 0.05\n")
        f.write("- **Multiple testing correction**: Benjamini-Hochberg FDR\n")
        f.write("- **Heatmap normalization**: Z-score by gene\n")
        f.write("- **PCA scaling**: Standard scaling (mean=0, variance=1)\n\n")
        
        f.write("## Comparison with Original Paper\n\n")
        f.write("The original paper (Li et al., 2025) reported:\n")
        f.write("- 223 shared DEGs between psoriasis and Crohn's disease\n")
        f.write("- Strong enrichment in cell cycle and immune response pathways\n")
        f.write("- Five key hub genes: KIF4A, DLGAP5, NCAPG, CCNB1, CEP55\n\n")
        
        f.write("Generated: " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
    
    print(f"  ✓ Created README file: {readme_path}")

if __name__ == "__main__":
    main()