import pandas as pd

# ==========================================
# STEP 1: Load the Top 500 Spreadsheet
# ==========================================
INPUT_CSV = "top500_7B7D_average.csv"  # Use your generated CSV from earlier
OUTPUT_CSV = "statistics_top_500_summary.csv"

# Paths to the full 310,000-molecule score files
FILE_FAR  = "../HTVS_runs/array_7B7D_far_all_mols/mega_sorted_docking_scores.txt"
FILE_HEAT = "../HTVS_runs/array_7B7D_HEAT_all_mols/mega_sorted_docking_scores.txt"
FILE_2IW3 = "../HTVS_runs/array_2IW3_all_mols/mega_sorted_docking_scores.txt"

print(f"Loading top 500 docking scores from {INPUT_CSV}...")
df_top500 = pd.read_csv(INPUT_CSV)

# ==========================================
# STEP 2: Define a Helper Function for Pairwise Math
# ==========================================

df_far_all  = pd.read_csv(FILE_FAR,  sep=r'\s+', header=None, names=['ZINC_ID', 'Score_7B7D_far'])
df_heat_all = pd.read_csv(FILE_HEAT, sep=r'\s+', header=None, names=['ZINC_ID', 'Score_7B7D_HEAT'])
df_2iw3_all = pd.read_csv(FILE_2IW3, sep=r'\s+', header=None, names=['ZINC_ID', 'Score_2IW3'])

# Merge all three full-library tables on ZINC_ID
df_full_library = df_far_all.merge(df_heat_all, on='ZINC_ID', how='inner') \
                            .merge(df_2iw3_all, on='ZINC_ID', how='inner')


def calculate_pocket_metrics(df_top, df_full, col_A, col_B, label_A, label_B):
    #Calculates MAD, MSD, Pearson r, and % agreement between two score columns.
    # 1. Vectorized subtraction across all 500 rows at once
    raw_diffs = df_top[col_A] - df_top[col_B]
    abs_diffs = raw_diffs.abs()
    
    # 2. Calculate Mean Absolute Difference (MAD)
    mad = abs_diffs.mean()
    
    # 3. Calculate Mean Signed Difference (MSD)
    msd = raw_diffs.mean()
    
    # 4. Calculate percentage of molecules within <= 1.0 kcal/mol
    # (abs_diffs <= 1.0) creates True/False (1/0); .mean() gives the decimal proportion
    pct_within_1 = (abs_diffs <= 1.0).mean() * 100.0

    
    return {
        "Comparison": f"{label_A} vs. {label_B}",
        "MAD (kcal/mol)": round(mad, 3),
        "MSD (kcal/mol)": round(msd, 3),
        "% within <= 1.0 kcal/mol": round(pct_within_1, 1)
    }

# ==========================================
# STEP 3: Analyze All Three Pairwise Comparisons
# ==========================================
print("Calculating distance metrics across pocket pairs...\n")

results = [
    # Pair 1: The two active-state 7B7D pockets
    calculate_pocket_metrics(
        df_top500, df_full_library,
        "Score_7B7D_far", "Score_7B7D_HEAT",
        "7B7D_far", "7B7D_HEAT"
    ),
    # Pair 2: Active (far) vs. Collapsed (2IW3)
    calculate_pocket_metrics(
        df_top500, df_full_library,
        "Score_7B7D_far", "Score_2IW3",
        "7B7D_far", "2IW3"
    ),
    # Pair 3: Active (HEAT) vs. Collapsed (2IW3)
    calculate_pocket_metrics(
        df_top500, df_full_library,
        "Score_7B7D_HEAT", "Score_2IW3",
        "7B7D_HEAT", "2IW3"
    )
]

# ==========================================
# STEP 4: Convert to Table, Display, and Export
# ==========================================
summary_df = pd.DataFrame(results)

# Print a formatted text table to the terminal
print(summary_df.to_string(index=False))
print("-" * 75)

# Save as a CSV to easily copy-paste into Excel or slides
summary_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSummary table saved to: {OUTPUT_CSV}")