# Save as: scripts/python/create_gse13355_metadata.py

import os
import pandas as pd
from ftplib import FTP
import gzip
import shutil
import time

def download_geo_metadata(geo_id):
    """Download SOFT file from GEO and extract metadata."""
    print(f"Downloading {geo_id} metadata...")
    
    # GEO FTP server details
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd(f'/geo/series/{geo_id[:-3]}nnn/{geo_id}/soft/')
    
    # List files
    files = ftp.nlst()
    soft_file = f"{geo_id}_family.soft.gz"
    
    if soft_file not in files:
        print(f"Warning: {soft_file} not found. Trying alternatives...")
        # Try other naming patterns
        for f in files:
            if 'soft' in f.lower() and 'gz' in f:
                soft_file = f
                break
    
    # Download the file
    local_file = f"{geo_id}_family.soft.gz"
    with open(local_file, 'wb') as f:
        ftp.retrbinary(f'RETR {soft_file}', f.write)
    
    ftp.quit()
    
    # Extract the gzipped file
    extracted_file = f"{geo_id}_family.soft"
    with gzip.open(local_file, 'rb') as f_in:
        with open(extracted_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    print(f"Downloaded and extracted {extracted_file}")
    return extracted_file

def parse_soft_file(soft_file):
    """Parse SOFT file to extract sample metadata."""
    print("Parsing SOFT file...")
    
    with open(soft_file, 'r') as f:
        content = f.read()
    
    # Split into sections
    sections = content.split('^')
    
    samples = []
    current_sample = {}
    in_sample_section = False
    
    for line in content.split('\n'):
        line = line.strip()
        
        if line.startswith('^SAMPLE ='):
            if current_sample:
                samples.append(current_sample)
            current_sample = {'GSM_ID': line.split('=')[1].strip()}
            in_sample_section = True
        
        elif in_sample_section and line.startswith('!Sample_'):
            if '=' in line:
                key = line.split('=')[0].replace('!Sample_', '').strip()
                value = line.split('=')[1].strip()
                current_sample[key] = value
            elif ':' in line:
                key = line.split(':')[0].replace('!Sample_', '').strip()
                value = line.split(':')[1].strip()
                current_sample[key] = value
        
        elif line.startswith('^') and not line.startswith('^SAMPLE'):
            if current_sample:
                samples.append(current_sample)
                current_sample = {}
            in_sample_section = False
    
    # Add the last sample
    if current_sample:
        samples.append(current_sample)
    
    df = pd.DataFrame(samples)
    print(f"Found {len(df)} samples")
    
    return df

def create_sample_groups(df, geo_id):
    """Create sample groups CSV based on metadata."""
    print("\nAvailable columns in metadata:")
    print(df.columns.tolist())
    
    # Look for columns that might indicate psoriasis vs control
    potential_group_cols = [col for col in df.columns 
                          if any(keyword in col.lower() for keyword in 
                                ['disease', 'diagnosis', 'group', 'characteristics', 'type', 'state', 'title', 'description'])]
    
    print("\nPotential group columns:")
    for col in potential_group_cols:
        print(f"  - {col}")
    
    # Try to automatically identify groups
    group_df = pd.DataFrame({'GSM_ID': df['GSM_ID']})
    group_df['Group'] = 'Unknown'
    
    # Common patterns for GSE13355 (psoriasis dataset)
    for col in potential_group_cols:
        if 'title' in col.lower() or 'characteristics' in col.lower():
            unique_values = df[col].dropna().unique()
            print(f"\nUnique values in {col}:")
            for val in unique_values[:10]:  # Show first 10
                print(f"  - {val}")
            
            # Auto-detect based on common terms
            for idx, row in df.iterrows():
                val = str(row[col]).lower()
                if 'psoriasis' in val or 'lesional' in val:
                    group_df.loc[group_df['GSM_ID'] == row['GSM_ID'], 'Group'] = 'Psoriasis'
                elif 'normal' in val or 'control' in val or 'healthy' in val:
                    group_df.loc[group_df['GSM_ID'] == row['GSM_ID'], 'Group'] = 'Control'
                elif 'non-lesional' in val:
                    group_df.loc[group_df['GSM_ID'] == row['GSM_ID'], 'Group'] = 'Exclude'  # Non-lesional
    
    # Count groups
    print(f"\nAuto-detected groups:")
    print(group_df['Group'].value_counts())
    
    # Save to CSV
    output_file = f"{geo_id}_sample_groups.csv"
    group_df.to_csv(output_file, index=False)
    print(f"\nSaved sample groups to: {output_file}")
    
    # Also save full metadata for reference
    full_metadata_file = f"{geo_id}_full_metadata.csv"
    df.to_csv(full_metadata_file, index=False)
    print(f"Saved full metadata to: {full_metadata_file}")
    
    return group_df

def main():
    geo_id = "GSE13355"
    
    try:
        # Download metadata
        soft_file = download_geo_metadata(geo_id)
        
        # Parse metadata
        metadata_df = parse_soft_file(soft_file)
        
        # Create sample groups
        groups_df = create_sample_groups(metadata_df, geo_id)
        
        print(f"\n=== SUMMARY ===")
        print(f"Total samples: {len(groups_df)}")
        print(f"Psoriasis samples: {sum(groups_df['Group'] == 'Psoriasis')}")
        print(f"Control samples: {sum(groups_df['Group'] == 'Control')}")
        print(f"Excluded samples: {sum(groups_df['Group'] == 'Exclude')}")
        print(f"Unknown samples: {sum(groups_df['Group'] == 'Unknown')}")
        
        print(f"\nNote: You may need to manually verify some samples.")
        print(f"Check the CSV file and update groups as needed.")
        
    except Exception as e:
        print(f"Error: {e}")
        print("\nAlternative: Creating template CSV based on paper information...")
        create_template_csv(geo_id)

def create_template_csv(geo_id):
    """Create a template CSV when automatic download fails."""
    print(f"\nCreating template CSV for {geo_id}...")
    
    # Based on paper: GSE13355 has 58 psoriasis and 64 control samples
    # We need to list all GSM IDs from the CEL files
    cel_path = f"../data/raw/{geo_id}"
    
    if os.path.exists(cel_path):
        cel_files = [f for f in os.listdir(cel_path) if f.endswith('.CEL.gz')]
        gsm_ids = [f.replace('.CEL.gz', '') for f in cel_files]
        
        df = pd.DataFrame({'GSM_ID': gsm_ids, 'Group': 'Unknown'})
        output_file = f"{geo_id}_sample_groups_template.csv"
        df.to_csv(output_file, index=False)
        
        print(f"Created template with {len(df)} samples")
        print(f"Saved to: {output_file}")
        print("\nYou need to manually annotate the groups based on:")
        print("1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13355")
        print("2. Check 'Sample title' column")
        print("3. Update CSV with 'Psoriasis' or 'Control'")
    else:
        print(f"Error: CEL file directory not found at {cel_path}")

if __name__ == "__main__":
    main()