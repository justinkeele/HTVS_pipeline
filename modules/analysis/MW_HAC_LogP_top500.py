import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# ==========================================
# STEP 0: Define File Paths
# ==========================================
TOP500_CSV       = "top500_7B7D_average.csv"
EXPERIMENTAL_SMI = "experimental_smiles.smi"
OUTPUT_CSV       = "MW_HAC_LogP_summary.csv"

# ==========================================
# STEP 1: Helper Function to Calculate RDKit Properties
# ==========================================
def compute_properties_from_smiles(smiles_list, id_list):
    """
    Takes a list of SMILES strings and returns a DataFrame of MW, HAC, and LogP.
    """
    records = []
    
    for smi, mol_id in zip(smiles_list, id_list):
        mol = Chem.MolFromSmiles(str(smi).strip())
        
        # Safety check: ensure RDKit successfully parsed the chemical structure
        if mol is not None:
            mw   = Descriptors.MolWt(mol)          # Average Molecular Weight (g/mol)
            hac  = mol.GetNumHeavyAtoms()          # Heavy Atom Count (non-hydrogen atoms)
            logp = Descriptors.MolLogP(mol)        # Wildman-Crippen LogP (lipophilicity)
            
            records.append({
                "ID": mol_id,
                "MW": mw,
                "HAC": hac,
                "LogP": logp
            })
        else:
            print(f"  [Warning] Could not parse SMILES for {mol_id}: {smi}")
            
    return pd.DataFrame(records)

# ==========================================
# STEP 2: Load and Process Top 500 Vina Hits
# ==========================================
print(f"Loading Top 500 computational hits from {TOP500_CSV}...")
df_top = pd.read_csv(TOP500_CSV)

# Extract SMILES and ZINC_ID from your existing spreadsheet
top_props_df = compute_properties_from_smiles(df_top['SMILES'], df_top['ZINC_ID'])

# ==========================================
# STEP 3: Load and Process Experimental Hits
# ==========================================
print(f"Loading Experimental hits from {EXPERIMENTAL_SMI}...")

# Read white-space delimited SMILES file (Col 0 = SMILES, Col 1 = Molecule Name/ID)
df_exp = pd.read_csv(EXPERIMENTAL_SMI, sep=r'\s+', header=None, names=['SMILES', 'EXP_ID'])
exp_props_df = compute_properties_from_smiles(df_exp['SMILES'], df_exp['EXP_ID'])

# ==========================================
# STEP 4: Build Summary Statistics Table
# ==========================================
def get_summary_row(prop_name, exp_series, top_series):
    """
    Calculates Mean ± SD, Median, Range, and Shift between Experimental and Top 500 sets.
    """
    exp_mean, exp_std = exp_series.mean(), exp_series.std()
    top_mean, top_std = top_series.mean(), top_series.std()
    
    delta_mean = top_mean - exp_mean
    
    return {
        "Property": prop_name,
        "Experimental (Mean ± SD)": f"{exp_mean:.1f} ± {exp_std:.1f}",
        "Experimental Median": f"{exp_series.median():.1f}",
        "Experimental Range": f"[{exp_series.min():.1f} – {exp_series.max():.1f}]",
        "Top 500 (Mean ± SD)": f"{top_mean:.1f} ± {top_std:.1f}",
        "Top 500 Median": f"{top_series.median():.1f}",
        "Top 500 Range": f"[{top_series.min():.1f} – {top_series.max():.1f}]",
        "Shift (Δ Mean)": f"{delta_mean:+.2f}"
    }

print("\nCalculating comparative physicochemical statistics...")

summary_rows = [
    get_summary_row("Molecular Weight (g/mol)", exp_props_df["MW"],   top_props_df["MW"]),
    get_summary_row("Heavy Atom Count (HAC)",   exp_props_df["HAC"],  top_props_df["HAC"]),
    get_summary_row("LogP (Lipophilicity)",     exp_props_df["LogP"], top_props_df["LogP"])
]

summary_df = pd.DataFrame(summary_rows)

# ==========================================
# STEP 5: Display and Export Results
# ==========================================
print("\n" + "="*90)
print(summary_df.to_string(index=False))
print("="*90)

summary_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved presentation-ready table to: {OUTPUT_CSV}")