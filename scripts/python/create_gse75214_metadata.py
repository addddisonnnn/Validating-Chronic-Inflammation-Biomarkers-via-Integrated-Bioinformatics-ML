import os
import pandas as pd

data_dir = "/projectnb/bf528/students/addisony/project1/Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper/data/raw/GSE75214"

files = os.listdir(data_dir)
cel_files = [f for f in files if f.endswith('.CEL.gz')]

metadata = []

for f in cel_files:
    gsm_id = f.split('_')[0]  # Extract GSM ID
    f_lower = f.lower()
    
    # Check if it's an ileum sample
    if 'ileum' in f_lower:
        if 'cd_ileum' in f_lower:
            group = 'CD'
        elif 'controle_ileum' in f_lower:
            group = 'Control'
        else:
            group = 'Exclude'  # Exclude UC ileum if any
        
        if group != 'Exclude':
            metadata.append({
                'GSM_ID': f.replace('.CEL.gz', ''),  # Use full filename without .CEL.gz
                'Group': group,
                'Tissue': 'ileum'
            })
    else:
        # Skip colon samples as per paper
        continue

df = pd.DataFrame(metadata)

# Count samples
print(f"Total ileum samples found: {len(df)}")
print(f"CD ileum samples: {sum(df['Group'] == 'CD')}")
print(f"Control ileum samples: {sum(df['Group'] == 'Control')}")

# Save to CSV
df[['GSM_ID', 'Group']].to_csv('GSE75214_sample_groups_fixed.csv', index=False)
print("\nCreated GSE75214_sample_groups_fixed.csv")
print("\nFirst few rows:")
print(df[['GSM_ID', 'Group']].head(10))