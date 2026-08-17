import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
    # Drop ATP and ADP, and keep only the HEAT pocket
    master_df = master_df[~master_df['Ligand'].isin(['ATP', 'ADP'])]
    master_df = master_df[master_df['Pocket'] == '7B7D_HEAT']

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
        
        # Calculate Consensus Score by averaging the three normalized raw scores
        group['Consensus'] = group[['Norm_CNN_Aff', 'Norm_CNN_Pose', 'Norm_Vina']].mean(axis=1)

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
                'X_Label': x_label,
                'Metric': metric_name,
                'BEDROC': bedroc_score,
                'NDCG': ndcg_score
            })

    return pd.DataFrame(results)

def plot_heatmaps(results_df):
    print("Generating Heatmap Grid...")
    
    # Set up a 2x2 grid
    sns.set_theme(style="white")
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True, sharey=True)
    
    # Define the 4 panels
    panels = [
        {'test_ver': '2', 'score_type': 'BEDROC', 'ax': axes[0, 0], 'title': 'Test 2 - BEDROC (Alpha=8)'},
        {'test_ver': '2', 'score_type': 'NDCG', 'ax': axes[0, 1], 'title': 'Test 2 - NDCG'},
        {'test_ver': '3', 'score_type': 'BEDROC', 'ax': axes[1, 0], 'title': 'Test 3 - BEDROC (Alpha=8)'},
        {'test_ver': '3', 'score_type': 'NDCG', 'ax': axes[1, 1], 'title': 'Test 3 - NDCG'}
    ]

    for p in panels:
        # Filter the data for this specific panel
        subset = results_df[results_df['Test_Version'] == p['test_ver']]
        
        # Pivot the data so X_Labels are columns and Metrics are rows
        matrix = subset.pivot(index='Metric', columns='X_Label', values=p['score_type'])
        
        # Reorder the matrix to match our desired layout perfectly
        matrix = matrix.reindex(index=Y_AXIS_ORDER, columns=X_AXIS_ORDER)

        # Draw the Heatmap. 
        # cmap="Blues" means Dark Blue = Good Score (1.0), Light Blue = Bad Score (0.0)
        sns.heatmap(
            matrix, ax=p['ax'], cmap="Blues", annot=True, fmt=".2f", 
            vmin=0.0, vmax=1.0, linewidths=1, linecolor='gray',
            cbar_kws={'label': p['score_type']}
        )
        
        p['ax'].set_title(p['title'], fontsize=14, weight='bold', pad=15)
        p['ax'].set_ylabel("")
        p['ax'].set_xlabel("")
        p['ax'].tick_params(axis='x', rotation=45)
        
        # THE VERTICAL LINE TRICK: Draw a thick white line to separate rescore from refinement
        # Because we added rescore_ex8, rescore now has 5 columns. We draw the line at index 5.
        p['ax'].axvline(x=5, color='black', linewidth=6)

    plt.tight_layout()
    output_img = os.path.join(OUTPUT_DIR, "HEAT_Heatmap_show_best_sorting_accuracy_test2_vs_3.svg")
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Success! Graphic saved to {output_img}")

if __name__ == '__main__':
    df_results = process_data()
    if not df_results.empty:
        plot_heatmaps(df_results)
    else:
        print("Error: No data successfully processed.")