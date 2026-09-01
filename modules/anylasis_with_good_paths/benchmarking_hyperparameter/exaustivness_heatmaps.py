import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from rdkit.ML.Scoring import Scoring

# 1. Configuration & Paths
INPUT_DIR = "../../../docs/graphs_with_good_paths"
OUTPUT_DIR = "../../../docs/graphs_with_good_paths"

# Define the exact layout we want for the X-axis and Y-axis
# ADDED: rescore_ex8 to the beginning of the list
X_AXIS_ORDER = [
    'rescore_ex8', 'rescore_ex16', 'rescore_ex32', 'rescore_ex64', 'rescore_ex128', 
    'refinement_ex8', 'refinement_ex16', 'refinement_ex32'
]
Y_AXIS_ORDER = ['Consensus', 'CNN_Affinity', 'CNN_Pose_Score', 'Vina_Affinity']

def assign_ground_truth(ligand_name):
    '''
    Assigns the binary active/inactive flag for BEDROC, 
    and the 0-3 weight scale for NDCG.
    '''
    if "tophit_high" in ligand_name:
        return 1, 3  # (Active for BEDROC, Weight 3 for NDCG)
    elif "tophit_low" in ligand_name:
        return 1, 2  # (Active for BEDROC, Weight 2 for NDCG)
    elif "no_inhibiton" in ligand_name:
        return 0, -0.5  # (Inactive for BEDROC, Weight -0.5 for NDCG)
    else:
        return 0, 0  # (Inactive for BEDROC, Weight 0 for Z-series junk)

def calculate_ndcg(true_weights):
    '''
    Calculates Normalized Discounted Cumulative Gain.
    Formula: DCG / Ideal_DCG
    '''
    # DCG rewards weights, but divides by log2 of the rank to penalize pushing hits down
    dcg = sum((2**w - 1) / np.log2(i + 2) for i, w in enumerate(true_weights))
    
    # Ideal DCG is what the score would be if the list was sorted perfectly
    ideal_weights = sorted(true_weights, reverse=True)
    idcg = sum((2**w - 1) / np.log2(i + 2) for i, w in enumerate(ideal_weights))
    
    return dcg / idcg if idcg > 0 else 0.0

# NEW FUNCTION: Min-Max Normalization
def min_max_normalize(series):
    '''
    Takes a column of scores and mathematically squashes them between 0.0 and 1.0.
    The worst score becomes exactly 0.0, and the best score becomes exactly 1.0.
    '''
    return (series - series.min()) / (series.max() - series.min())

def process_data():
    # 1. Load the Data
    print("Loading master CSVs...")
    df_cnn_aff = pd.read_csv(os.path.join(INPUT_DIR, "master_cnn_affinity.csv"))
    df_cnn_pose = pd.read_csv(os.path.join(INPUT_DIR, "master_cnn_pose_score.csv"))
    df_vina = pd.read_csv(os.path.join(INPUT_DIR, "master_vina_affinity.csv"))

    # Rename the 'Score' column in each so we can merge them safely
    df_cnn_aff = df_cnn_aff.rename(columns={'Score': 'CNN_Affinity'})
    df_cnn_pose = df_cnn_pose.rename(columns={'Score': 'CNN_Pose_Score'})
    df_vina = df_vina.rename(columns={'Score': 'Vina_Affinity'})

    # Merge them into one big table
    master_df = df_cnn_aff.merge(df_cnn_pose, on=['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Test_Version'])
    master_df = master_df.merge(df_vina, on=['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Test_Version'])

    # 2. Filter & Clean
    # Drop ATP and ADP; keep all pockets so we can plot multiple pockets as rows
    master_df = master_df[~master_df['Ligand'].isin(['ATP', 'ADP'])]

    # NEW: Standardize Vina Affinity direction by multiplying by -1
    # Now, higher is better for ALL metrics.
    master_df['Vina_Affinity'] = master_df['Vina_Affinity'] * -1

    # Apply the Ground Truth weights
    master_df[['Is_Active', 'NDCG_Weight']] = master_df['Ligand'].apply(assign_ground_truth).apply(pd.Series)

    results = []

    # 3. Process each unique run mathematically
    # Group by the specific run parameters
    for name, group in master_df.groupby(['Run_Name', 'Mode', 'Exhaustiveness', 'Test_Version']):
        run_name, mode, ex, test_ver = name
        x_label = f"{mode}_{ex}"
        
        # Skip runs that aren't in our planned X-axis (like CPU runs or weird exhaustiveness)
        if x_label not in X_AXIS_ORDER:
            continue

        group = group.copy()
        
        # NEW: Normalize the raw scores to a 0-1 scale so they can be averaged fairly
        group['Norm_CNN_Aff'] = min_max_normalize(group['CNN_Affinity'])
        group['Norm_CNN_Pose'] = min_max_normalize(group['CNN_Pose_Score'])
        group['Norm_Vina'] = min_max_normalize(group['Vina_Affinity'])
        
        # Calculate Weighted Consensus Score:
        # Give more weight to Pose (0.45), moderate weight to Vina (0.30), and less to Affinity (0.25)
        group['Consensus'] = (
            (0.25 * group['Norm_CNN_Aff']) + 
            (0.45 * group['Norm_CNN_Pose']) + 
            (0.30 * group['Norm_Vina'])
        )

        # 4. Calculate BEDROC and NDCG for each metric
        # Because we flipped Vina, ALL metrics now sort Descending (False) where Higher = Better
        metrics = [
            'Consensus', 
            'CNN_Affinity', 
            'CNN_Pose_Score', 
            'Vina_Affinity'
        ]

        for metric_name in metrics:
            sorted_group = group.sort_values(by=metric_name, ascending=False)
            
            # --- NDCG Math ---
            true_weights = sorted_group['NDCG_Weight'].tolist()
            ndcg_score = calculate_ndcg(true_weights)
            
            # --- RDKit BEDROC Math ---
            # RDKit CalcBEDROC requires a list of lists: [[active_flag], [active_flag], ...]
            # It assumes the list is ALREADY sorted from best score to worst score.
            rdkit_data = [[flag] for flag in sorted_group['Is_Active'].tolist()]
            bedroc_score = Scoring.CalcBEDROC(rdkit_data, 0, alpha=8.0)

            results.append({
                'Test_Version': str(test_ver),
                'Pocket': group['Pocket'].iloc[0],
                'X_Label': x_label,
                'Metric': metric_name,
                'BEDROC': bedroc_score,
                'NDCG': ndcg_score
            })

    return pd.DataFrame(results)

def plot_heatmaps(results_df):
    print("Generating Heatmap Grid...")
    
    # Set up a 3x2 grid: rows are pockets (HEAT, far, 2IW3), cols are [BEDROC, NDCG]
    sns.set_theme(style="white")
    pockets = [
        ('7B7D_HEAT', 'HEAT'),
        ('7B7D_far', 'Far'),
        ('2IW3', '2IW3')
    ]
    score_types = ['BEDROC', 'NDCG']

    fig, axes = plt.subplots(len(pockets), len(score_types), figsize=(16, 12), sharex=True, sharey=True)

    # --- THE FIX: BUILD THE UNIFIED COLORMAP ONCE ---
    blues = plt.get_cmap('Blues')
    blue_samples = [blues(t) for t in np.linspace(0.0, 1.0, 128)]
    dark_red = np.array([0.55, 0.0, 0.0])
    white = np.array([1.0, 1.0, 1.0])
    red_samples = [tuple(dark_red * (1.0 - t) + white * t) for t in np.linspace(0.0, 1.0, 128)]
    
    color_list = red_samples + [tuple(white)] + blue_samples
    diverging_cmap = LinearSegmentedColormap.from_list('custom_diverging', color_list)
    # ------------------------------------------------

    # Pre-select rows: for each Pocket/X_Label/Metric pick the row with the highest Test_Version
    df = results_df.copy()
    df['Test_Version_int'] = df['Test_Version'].astype(int)

    # Exclude Test_Version 1 rows for HEAT and far pockets so those cells remain blank
    df_for_sel = df[~((df['Pocket'].isin(['7B7D_HEAT', '7B7D_far'])) & (df['Test_Version_int'] == 1))]

    # Pick the latest Test_Version per Pocket/X_Label/Metric from the filtered set
    if not df_for_sel.empty:
        idx = df_for_sel.groupby(['Pocket', 'X_Label', 'Metric'])['Test_Version_int'].idxmax()
        df_sel = df.loc[idx]
    else:
        df_sel = pd.DataFrame(columns=df.columns)

    for i, (pocket_key, pocket_label) in enumerate(pockets):
        for j, score_type in enumerate(score_types):
            ax = axes[i, j]

            # Filter the selected data for this pocket
            subset = df_sel[df_sel['Pocket'] == pocket_key]

            # Pivot the data so X_Labels are columns and Metrics are rows
            matrix = subset.pivot(index='Metric', columns='X_Label', values=score_type)

            # Reorder the matrix to match our desired layout; keep NaNs as-is
            matrix = matrix.reindex(index=Y_AXIS_ORDER, columns=X_AXIS_ORDER)

            # --- THE FIX: APPLY TWO-SLOPE NORMALIZATION ---
            if score_type == 'BEDROC':
                # BEDROC: bounds are 0.0 to 1.0, random chance is 0.42
                norm = TwoSlopeNorm(vmin=0.0, vcenter=0.42, vmax=1.0)
            else:
                # NDCG: bounds are 0.424 to 1.0, random chance is 0.644
                norm = TwoSlopeNorm(vmin=0.424, vcenter=0.644, vmax=1.0)
                
            sns.heatmap(
                matrix, ax=ax, cmap=diverging_cmap, norm=norm, annot=True, fmt=".2f",
                linewidths=1, linecolor='gray', cbar_kws={'label': score_type}
            )
            # ----------------------------------------------

            ax.set_title(f"{pocket_label} - {score_type}", fontsize=12, weight='bold', pad=12)
            ax.set_ylabel("")
            ax.set_xlabel("")
            ax.tick_params(axis='x', rotation=45)

            # Vertical separator between rescore and refinement columns (after 5 rescore columns)
            ax.axvline(x=5, color='black', linewidth=6)

    plt.tight_layout()
    output_img = os.path.join(OUTPUT_DIR, "HEAT_Heatmap_show_best_sorting_accuracy_all_three_pockets_weighted_consensus2.png")
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Success! Graphic saved to {output_img}")

if __name__ == '__main__':
    df_results = process_data()
    if not df_results.empty:
        plot_heatmaps(df_results)
    else:
        print("Error: No data successfully processed.")