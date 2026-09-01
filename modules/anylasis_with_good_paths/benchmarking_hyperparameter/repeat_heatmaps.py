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

# Define the exact 9-column layout for the X-axis across the three pockets
X_AXIS_ORDER = [
    'HEAT_ex32', 'HEAT_ex64', 'HEAT_ex128',
    'Far_ex32', 'Far_ex64', 'Far_ex128',
    '2IW3_ex32', '2IW3_ex64', '2IW3_ex128'
]

Y_AXIS_ORDER = ['Consensus', 'CNN_Affinity', 'CNN_Pose_Score', 'Vina_Affinity']

def assign_ground_truth(ligand_name):
    if "tophit_high" in ligand_name: return 1, 3
    elif "tophit_low" in ligand_name: return 1, 2
    elif "no_inhibiton" in ligand_name: return 0, -0.5
    else: return 0, 0

def calculate_ndcg(true_weights):
    dcg = sum((2**w - 1) / np.log2(i + 2) for i, w in enumerate(true_weights))
    ideal_weights = sorted(true_weights, reverse=True)
    idcg = sum((2**w - 1) / np.log2(i + 2) for i, w in enumerate(ideal_weights))
    return dcg / idcg if idcg > 0 else 0.0

def min_max_normalize(series):
    if series.max() == series.min(): return series * 0.0 + 0.5 
    return (series - series.min()) / (series.max() - series.min())

def process_data():
    print("Loading 10-replicate CSVs...")
    df_cnn_aff = pd.read_csv(os.path.join(INPUT_DIR, "repeat_runs_cnn_affinity.csv")).rename(columns={'Score': 'CNN_Affinity'})
    df_cnn_pose = pd.read_csv(os.path.join(INPUT_DIR, "repeat_runs_cnn_pose_score.csv")).rename(columns={'Score': 'CNN_Pose_Score'})
    df_vina = pd.read_csv(os.path.join(INPUT_DIR, "repeat_runs_vina_affinity.csv")).rename(columns={'Score': 'Vina_Affinity'})

    # Merge on the new Repeat_Number column
    merge_cols = ['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Repeat_Number']
    master_df = df_cnn_aff.merge(df_cnn_pose, on=merge_cols).merge(df_vina, on=merge_cols)

    # Clean the data
    master_df = master_df[~master_df['Ligand'].isin(['ATP', 'ADP'])]
    master_df['Vina_Affinity'] = master_df['Vina_Affinity'] * -1
    master_df[['Is_Active', 'NDCG_Weight']] = master_df['Ligand'].apply(assign_ground_truth).apply(pd.Series)

    # Create shortened pocket labels to keep the graph axis clean
    pocket_map = {'7B7D_HEAT': 'HEAT', '7B7D_far': 'Far', '2IW3': '2IW3'}
    master_df['Short_Pocket'] = master_df['Pocket'].map(pocket_map)

    # --- THE FIX: ENSEMBLE AVERAGING ---
    # We group by the specific ligand and setup, and calculate the mean score across all 10 repeats!
    print("Averaging raw scores across 10 replicates...")
    agg_dict = {
        'CNN_Affinity': 'mean',
        'CNN_Pose_Score': 'mean',
        'Vina_Affinity': 'mean',
        'Is_Active': 'first',      # The ground truth doesn't change, just grab the first one
        'NDCG_Weight': 'first'
    }
    # This completely eliminates Repeat_Number, leaving 1 robust score per molecule
    avg_df = master_df.groupby(['Ligand', 'Short_Pocket', 'Exhaustiveness']).agg(agg_dict).reset_index()

    results = []

    # Now we loop through our clean, averaged configurations
    for name, group in avg_df.groupby(['Short_Pocket', 'Exhaustiveness']):
        short_pocket, ex = name

        # --- THE FIX ---
        # Removed the extra '_ex' string injection
        x_label = f"{short_pocket}_{ex}"
        
        if x_label not in X_AXIS_ORDER:
            continue

        group = group.copy()
        
        # Normalize the ENSEMBLE AVERAGED scores
        group['Norm_CNN_Aff'] = min_max_normalize(group['CNN_Affinity'])
        group['Norm_CNN_Pose'] = min_max_normalize(group['CNN_Pose_Score'])
        group['Norm_Vina'] = min_max_normalize(group['Vina_Affinity'])
        
        # Weighted Consensus: 25% Aff, 45% Pose, 30% Vina
        group['Consensus'] = (
            (0.25 * group['Norm_CNN_Aff']) + 
            (0.45 * group['Norm_CNN_Pose']) + 
            (0.30 * group['Norm_Vina'])
        )

        metrics = ['Consensus', 'CNN_Affinity', 'CNN_Pose_Score', 'Vina_Affinity']

        for metric_name in metrics:
            sorted_group = group.sort_values(by=metric_name, ascending=False)
            
            ndcg_score = calculate_ndcg(sorted_group['NDCG_Weight'].tolist())
            rdkit_data = [[flag] for flag in sorted_group['Is_Active'].tolist()]
            bedroc_score = Scoring.CalcBEDROC(rdkit_data, 0, alpha=8.0)

            results.append({
                'X_Label': x_label,
                'Metric': metric_name,
                'BEDROC': bedroc_score,
                'NDCG': ndcg_score
            })

    return pd.DataFrame(results)

def plot_master_heatmap(results_df):
    print("Generating Master Averaged Heatmap...")
    
    sns.set_theme(style="white")
    
    # 2 rows, 1 column. sharex=True locks their X-axes perfectly together.
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    score_types = ['BEDROC', 'NDCG']

    # --- THE FIX: BUILD THE UNIFIED COLORMAP ONCE ---
    blues = plt.get_cmap('Blues')
    blue_samples = [blues(t) for t in np.linspace(0.0, 1.0, 128)]
    dark_red = np.array([0.55, 0.0, 0.0])
    white = np.array([1.0, 1.0, 1.0])
    red_samples = [tuple(dark_red * (1.0 - t) + white * t) for t in np.linspace(0.0, 1.0, 128)]
    
    color_list = red_samples + [tuple(white)] + blue_samples
    diverging_cmap = LinearSegmentedColormap.from_list('custom_diverging', color_list)
    # ------------------------------------------------

    for j, score_type in enumerate(score_types):
        ax = axes[j]

        # Pivot the data so X_Labels are columns and Metrics are rows
        matrix = results_df.pivot(index='Metric', columns='X_Label', values=score_type)

        # Reorder the matrix perfectly
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

        ax.set_title(f"10-Run Average Sorting Accuracy - {score_type}", fontsize=16, weight='bold', pad=15)
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.tick_params(axis='x', rotation=45, labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        # --- THE VISUAL DELIMITERS ---
        # Draw thick black lines to separate HEAT, Far, and 2IW3
        ax.axvline(x=3, color='black', linewidth=5)
        ax.axvline(x=6, color='black', linewidth=5)

    plt.tight_layout()
    output_img = os.path.join(OUTPUT_DIR, "Master_Ensemble_Heatmap2.png")
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Success! Master Heatmap saved to {output_img}")

if __name__ == '__main__':
    df_results = process_data()
    if not df_results.empty:
        plot_master_heatmap(df_results)
    else:
        print("Error: No data successfully processed.")