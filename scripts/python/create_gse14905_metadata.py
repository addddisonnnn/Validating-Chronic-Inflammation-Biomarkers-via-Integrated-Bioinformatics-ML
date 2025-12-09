import os
import pandas as pd

data_dir = "/projectnb/bf528/students/addisony/project1/Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper/data/raw/GSE14905"

files = os.listdir(data_dir)
cel_files = [f for f in files if f.endswith('.CEL.gz')]

metadata = []

# Based on GEO annotations:
normal_samples = list(range(372286, 372307))  # GSM372286 to GSM372306
lesional_samples = [
    372308, 372310, 372312, 372314, 372316, 372318, 372320, 372322, 372324, 372325,
    372327, 372329, 372331, 372333, 372334, 372336, 372338, 372340, 372342, 372344,
    372346, 372348, 372350, 372352, 372354, 372356, 372358, 372360, 372362, 372364,
    372365, 372366, 372367
]

for f in cel_files:
    gsm_id = f.replace('.CEL.gz', '')
    gsm_num = int(gsm_id.replace('GSM', ''))
    
    if gsm_num in normal_samples:
        group = 'Control'
    elif gsm_num in lesional_samples:
        group = 'Psoriasis'
    else:
        group = 'Exclude'  # Uninvolved/non-lesional samples
    
    if group != 'Exclude':
        metadata.append({
            'GSM_ID': gsm_id,
            'Group': group
        })

df = pd.DataFrame(metadata)

# Count samples
print(f"Total samples to include: {len(df)}")
print(f"Control samples: {sum(df['Group'] == 'Control')}")
print(f"Psoriasis Lesional samples: {sum(df['Group'] == 'Psoriasis')}")

# Save to CSV
df.to_csv('GSE14905_sample_groups.csv', index=False)
print("\nCreated GSE14905_sample_groups.csv")
print("\nFirst few rows:")
print(df.head(10))