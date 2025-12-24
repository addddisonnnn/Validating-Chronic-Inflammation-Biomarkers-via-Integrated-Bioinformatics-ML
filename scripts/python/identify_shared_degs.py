# Save as: scripts/python/identify_shared_degs_fixed.py

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os

# Set the base directory
base_dir = "/projectnb/bf528/students/addisony/project1/Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper"

def load_degs():
    """Load DEG files from the processed directory."""
    deg_files = {
        'psoriasis_train': os.path.join(base_dir, 'data/processed/GSE13355_DEGs.csv'),
        'crohns_train': os.path.join(base_dir, 'data/processed/GSE75214_DEGs.csv'),
        'psoriasis_val': os.path.join(base_dir, 'data/processed/GSE14905_DEGs.csv'),
        'crohns_val': os.path.join(base_dir, 'data/processed/GSE102133_DEGs.csv')
    }
    
    deg_data = {}
    for name, path in deg_files.items():
        print(f"Loading {name} from {path}")
        
        if os.path.exists(path):
            df = pd.read_csv(path)
            print(f"  File loaded: {df.shape[0]} rows, {df.shape[1]} columns")
            print(f"  Columns: {df.columns.tolist()}")
            
            # Find the gene column
            gene_col = None
            possible_names = ['Gene', 'gene', 'SYMBOL', 'symbol', 'gene_symbol']
            
            for col in df.columns:
                if any(possible in str(col).lower() for possible in ['gene', 'symbol']):
                    gene_col = col
                    break
            
            if gene_col:
                genes = set(df[gene_col].dropna().astype(str).str.strip())
                deg_data[name] = genes
                print(f"  Found {len(genes)} genes using column '{gene_col}'")
            else:
                # Try using first column if it looks like gene names
                first_col = df.columns[0]
                genes = set(df[first_col].dropna().astype(str).str.strip())
                deg_data[name] = genes
                print(f"  Using first column '{first_col}': {len(genes)} genes")
        else:
            print(f"  ERROR: File not found at {path}")
            deg_data[name] = set()
    
    return deg_data

def find_shared_degs(deg_data):
    """Find shared DEGs between datasets."""
    print("\n" + "="*60)
    print("FINDING SHARED DEGs")
    print("="*60)
    
    # Shared between training datasets
    psoriasis_train = deg_data.get('psoriasis_train', set())
    crohns_train = deg_data.get('crohns_train', set())
    
    shared_train = psoriasis_train.intersection(crohns_train)
    print(f"\nPsoriasis (GSE13355): {len(psoriasis_train)} genes")
    print(f"Crohn's (GSE75214): {len(crohns_train)} genes")
    print(f"Shared (training sets): {len(shared_train)} genes")
    
    # Shared between validation datasets
    psoriasis_val = deg_data.get('psoriasis_val', set())
    crohns_val = deg_data.get('crohns_val', set())
    
    shared_val = psoriasis_val.intersection(crohns_val)
    print(f"\nPsoriasis validation (GSE14905): {len(psoriasis_val)} genes")
    print(f"Crohn's validation (GSE102133): {len(crohns_val)} genes")
    print(f"Shared (validation sets): {len(shared_val)} genes")
    
    # Shared across all
    shared_all = shared_train.intersection(shared_val)
    print(f"\nShared across all 4 datasets: {len(shared_all)} genes")
    
    return {
        'shared_train': shared_train,
        'shared_val': shared_val,
        'shared_all': shared_all,
        'psoriasis_train': psoriasis_train,
        'crohns_train': crohns_train
    }

def create_venn_diagram(psoriasis_genes, crohns_genes, output_path):
    """Create Venn diagram."""
    print(f"\nCreating Venn diagram...")
    
    if len(psoriasis_genes) == 0 or len(crohns_genes) == 0:
        print("ERROR: Cannot create Venn diagram - one or both gene sets are empty")
        return
    
    plt.figure(figsize=(10, 8))
    
    # Create the Venn diagram
    v = venn2([psoriasis_genes, crohns_genes], 
              set_labels=('Psoriasis\n(GSE13355)', 'Crohn\'s Disease\n(GSE75214)'))
    
    # Customize the diagram
    if v.get_patch_by_id('10'):
        v.get_patch_by_id('10').set_color('lightcoral')
    if v.get_patch_by_id('01'):
        v.get_patch_by_id('01').set_color('lightblue')
    if v.get_patch_by_id('11'):
        v.get_patch_by_id('11').set_color('lightgreen')
    
    # Update labels with counts
    v.get_label_by_id('10').set_text(f'Psoriasis only\n{len(psoriasis_genes - crohns_genes)}')
    v.get_label_by_id('01').set_text(f'CD only\n{len(crohns_genes - psoriasis_genes)}')
    v.get_label_by_id('11').set_text(f'Shared\n{len(psoriasis_genes & crohns_genes)}')
    
    plt.title('Shared DEGs: Psoriasis vs Crohn\'s Disease (Training Sets)', 
              fontsize=16, fontweight='bold', pad=20)
    
    # Save the figure
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Venn diagram saved to: {output_path}")

def save_shared_genes(shared_data, output_dir):
    """Save shared genes to CSV files."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Save training shared genes
    if shared_data['shared_train']:
        train_df = pd.DataFrame({'Gene': sorted(list(shared_data['shared_train']))})
        train_path = os.path.join(output_dir, 'shared_degs_training.csv')
        train_df.to_csv(train_path, index=False)
        print(f"Saved training shared genes: {train_path} ({len(train_df)} genes)")
    
    # Save validation shared genes
    if shared_data['shared_val']:
        val_df = pd.DataFrame({'Gene': sorted(list(shared_data['shared_val']))})
        val_path = os.path.join(output_dir, 'shared_degs_validation.csv')
        val_df.to_csv(val_path, index=False)
        print(f"Saved validation shared genes: {val_path} ({len(val_df)} genes)")
    
    # Save all shared genes
    if shared_data['shared_all']:
        all_df = pd.DataFrame({'Gene': sorted(list(shared_data['shared_all']))})
        all_path = os.path.join(output_dir, 'shared_degs_all.csv')
        all_df.to_csv(all_path, index=False)
        print(f"Saved all shared genes: {all_path} ({len(all_df)} genes)")

def create_summary_table(deg_data, shared_data):
    """Create summary statistics table."""
    summary_data = []
    
    # Training datasets
    summary_data.append({
        'Dataset': 'Psoriasis (GSE13355)',
        'Type': 'Training',
        'Total DEGs': len(deg_data['psoriasis_train']),
        'Shared DEGs': len(shared_data['shared_train']),
        'Percent Shared': f"{len(shared_data['shared_train'])/len(deg_data['psoriasis_train'])*100:.1f}%"
    })
    
    summary_data.append({
        'Dataset': 'Crohn\'s (GSE75214)',
        'Type': 'Training',
        'Total DEGs': len(deg_data['crohns_train']),
        'Shared DEGs': len(shared_data['shared_train']),
        'Percent Shared': f"{len(shared_data['shared_train'])/len(deg_data['crohns_train'])*100:.1f}%"
    })
    
    # Validation datasets
    summary_data.append({
        'Dataset': 'Psoriasis (GSE14905)',
        'Type': 'Validation',
        'Total DEGs': len(deg_data['psoriasis_val']),
        'Shared DEGs': len(shared_data['shared_val']),
        'Percent Shared': f"{len(shared_data['shared_val'])/len(deg_data['psoriasis_val'])*100:.1f}%"
    })
    
    summary_data.append({
        'Dataset': 'Crohn\'s (GSE102133)',
        'Type': 'Validation',
        'Total DEGs': len(deg_data['crohns_val']),
        'Shared DEGs': len(shared_data['shared_val']),
        'Percent Shared': f"{len(shared_data['shared_val'])/len(deg_data['crohns_val'])*100:.1f}%"
    })
    
    summary_df = pd.DataFrame(summary_data)
    summary_path = os.path.join(base_dir, 'results', 'shared_degs_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(summary_df.to_string(index=False))
    print(f"\nSummary saved to: {summary_path}")
    
    return summary_df

def main():
    print("="*60)
    print("IDENTIFYING SHARED DEGs: Psoriasis vs Crohn's Disease")
    print("="*60)
    
    # Load DEG results
    deg_data = load_degs()
    
    # Find shared DEGs
    shared_data = find_shared_degs(deg_data)
    
    # Create Venn diagram (if we have data)
    if shared_data['psoriasis_train'] and shared_data['crohns_train']:
        venn_path = os.path.join(base_dir, 'results', 'plots', 'shared_degs_venn.png')
        create_venn_diagram(shared_data['psoriasis_train'], 
                           shared_data['crohns_train'], 
                           venn_path)
    
    # Save shared genes
    results_dir = os.path.join(base_dir, 'results')
    save_shared_genes(shared_data, results_dir)
    
    # Create summary table
    summary_df = create_summary_table(deg_data, shared_data)
    
    # Compare with paper
    print("\n" + "="*60)
    print("COMPARISON WITH PAPER")
    print("="*60)
    print("Paper reported: 223 shared DEGs between GSE13355 and GSE75214")
    print(f"Our analysis: {len(shared_data['shared_train'])} shared DEGs")
    
    if len(shared_data['shared_train']) > 0:
        print("\nFirst 20 shared genes:")
        shared_list = sorted(list(shared_data['shared_train']))
        for i, gene in enumerate(shared_list[:20], 1):
            print(f"  {i:2d}. {gene}")

if __name__ == "__main__":
    main()