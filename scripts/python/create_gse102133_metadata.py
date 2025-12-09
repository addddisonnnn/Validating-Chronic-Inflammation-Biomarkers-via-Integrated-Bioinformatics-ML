import os
import pandas as pd

data_dir = "/projectnb/bf528/students/addisony/project1/Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper/data/raw/GSE102133"

files = os.listdir(data_dir)
cel_files = [f for f in files if f.endswith('.CEL.gz')]

print(f"Total CEL files: {len(cel_files)}")

metadata = []
for f in cel_files:
    gsm_id = f.replace('.CEL.gz', '')
    f_lower = f.lower()
    
    if 'control' in f_lower:
        group = 'Control'
    elif 'cd' in f_lower:
        group = 'CD'
    else:
        group = 'Unknown'
    
    metadata.append({
        'GSM_ID': gsm_id,
        'Filename': f,
        'Group': group
    })

df = pd.DataFrame(metadata)

# Count samples
print(f"\nSample counts:")
print(df['Group'].value_counts())

# Save to CSV
df[['GSM_ID', 'Group']].to_csv('GSE102133_sample_groups.csv', index=False)
print("\nCreated GSE102133_sample_groups.csv")
print("\nFirst few rows:")
print(df[['GSM_ID', 'Group']].head(10))
