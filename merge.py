import pandas as pd
import numpy as np

# Constants
R = 8.314 / 4184  # Gas constant 
T = 298.15  # Temperature in Kelvin

# Load UniProt-PDB mapping (TSV format)
mapping_df = pd.read_csv("idmapping_(uniprot-pdb).tsv", sep="\t", names=["UniProt_ID", "PDB"])

#mapping_df = pd.read_csv("uniprot_pdb.tsv", sep="\t", names=["UniProt_ID", "PDB"])

#print(mapping_df.head())

# Expand multiple PDB entries into separate rows
#mapping_df = mapping_df.assign(PDBs=mapping_df["PDB"].str.split(";")).explode("PDB")
#mapping_df["PDB"] = mapping_df["PDB"].str[:4].str.upper()

# Load SKEMPI dataset
skempi_df = pd.read_csv("skempi_v2.csv", delimiter=';')

#skempi_df.to_csv("skempi_try.csv", index=False)

# Remove rows where either 'affinity_wt' or 'affinity_mut' is missing
skempi_df = skempi_df.dropna(subset=['Affinity_mut (M)', 'Affinity_wt (M)'], how="any")

#print(skempi_df['Affinity_mut (M)'].unique())
#print(skempi_df['Affinity_wt (M)'].unique())

# Convert parsed affinity columns to numeric
skempi_df['Affinity_wt (M)'] = pd.to_numeric(skempi_df['Affinity_wt (M)'], errors='coerce')
skempi_df['Affinity_mut (M)'] = pd.to_numeric(skempi_df['Affinity_mut (M)'], errors='coerce')

epsilon = 1e-10

# Calculate ΔG for wild-type and mutant
skempi_df['deltaG_mut'] = R * T * np.log(skempi_df['Affinity_mut (M)'] + epsilon)
skempi_df['deltaG_wt'] = R * T * np.log(skempi_df['Affinity_wt (M)'] + epsilon)

# Calculate ΔΔG
skempi_df['delta_deltaG'] = skempi_df['deltaG_mut'] - skempi_df['deltaG_wt']

# Extract the 4-character PDB ID and convert to uppercase
skempi_df["#Pdb"] = skempi_df["#Pdb"].str[:4].str.upper()

# Function to parse mutation string
def parse_mutation(mutation):
    """Extracts wild-type amino acid, position, and mutant amino acid from mutation string."""
    try:
        # Extract wild-type amino acid (first two characters)
        wt = mutation[:2]
        
        # Extract position (digits in the middle)
        pos_str = ""
        for char in mutation[2:]:
            if char.isdigit():
                pos_str += char
            else:
                break
        pos = int(pos_str) if pos_str else None
        
        # Extract mutant amino acid (last character)
        mut = mutation[-1]
        
        return wt, pos, mut
    except Exception as e:
        print(f"Error parsing mutation: {mutation}. Error: {e}")
        return None, None, None

# Add parsed mutation columns to SKEMPI data
parsed_mutations = skempi_df["Mutation(s)_cleaned"].apply(parse_mutation)

# Check for parsing errors
if any(mutation == (None, None, None) for mutation in parsed_mutations):
    print("Warning: Some mutations could not be parsed. Check the mutation strings.")

# Convert parsed mutations to DataFrame
skempi_df[["wt", "pos", "mut"]] = pd.DataFrame(parsed_mutations.tolist(), index=skempi_df.index)

# Simplify mapping: Use the first character of the two-letter code
skempi_df["wt"] = skempi_df["wt"].str[0]

#print(skempi_df.head())

# Filter and rename SKEMPI columns
skempi_filtered = skempi_df[["#Pdb", "wt", "pos", "mut", "delta_deltaG"]].rename(
    columns={
        "#Pdb": "PDB_ID",
        "wt": "wt",
        "pos": "Mutation_Position",
        "mut": "mut",
        "delta_deltaG": "delta_deltaG" 
    }
)

#print(skempi_filtered.head())

skempi_filtered.to_csv("skempi_filtered.csv", index=False)

mapping_df2 = pd.read_csv("idmapping_(pdb-uniprot).tsv", sep="\t", usecols=[0, 1], names=["PDB", "Uniprot_ID"])

# Merge SKEMPI with UniProt mapping
skempi_uniprot_df = skempi_filtered.merge(mapping_df2, left_on="PDB_ID", right_on="PDB", how="inner")

#print(skempi_uniprot_df.head())

#skempi_uniprot_df.to_csv("skempi_final.csv", index=False)

# Load PPI dataset
ppi_df = pd.read_csv("ppi.csv")

#ppi_pdb_df = ppi_df.merge(mapping_df, left_on= "uniprot_id", right_on= "UniProt_ID", how= "inner")

# Merge SKEMPI-UniProt data with PPI dataset
ppi_ddg_df = ppi_df.merge(
    skempi_uniprot_df,
    left_on=["uniprot_id", "aa_ProtPosition", "sequence"],
    right_on=["Uniprot_ID", "Mutation_Position", "wt"],
    how="inner"
)

# Save the final dataset
ppi_ddg_df.to_csv("ppi_ddg_dataset.csv", index=False)

