import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.ML.Scoring import Scoring
import matplotlib.ticker as ticker

# 1. Configuration & Paths
INPUT_DIR = "../../../docs/graphs_with_good_paths"
OUTPUT_DIR = "../../../docs/graphs_with_good_paths"

# Define the exact order of the rows for our graph
POCKETS = ['7B7D_HEAT', '7B7D_far', '2IW3']

# --- HELPER FUNCTIONS ---
def clean_base_config(run_name):
    '''
    Strips away the _runN and _test_N suffixes so we can perfectly match 
    hardware logs to scoring logs based on their "Base Configuration".
    '''
    base = re.sub(r'_run\d+$', '', run_name)
    base = re.sub(r'_test_\d+$', '', base)
    base = re.sub(r'_\d+$', '', base)
    return base

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

# --- MAIN LOGIC ---
def process_data():
    print("Loading Telemetry & Scoring Data...")
    
    # 1. Process Hardware Data (Averaging Time)
    df_hw = pd.read_csv(os.path.join(INPUT_DIR, "hardware_telemetry_cleanup.csv"))
    df_hw['Base_Config'] = df_hw['Run_Name'].apply(clean_base_config)
    hw_avg = df_hw.groupby('Base_Config')['Avg_Time_Per_Mol_s'].mean().reset_index()
    hw_avg['Time_Min'] = hw_avg['Avg_Time_Per_Mol_s'] / 60.0

    # 2. Process Scoring Data
    df_cnn_aff = pd.read_csv(os.path.join(INPUT_DIR, "master_cnn_affinity.csv")).rename(columns={'Score': 'CNN_Affinity'})
    df_cnn_pose = pd.read_csv(os.path.join(INPUT_DIR, "master_cnn_pose_score.csv")).rename(columns={'Score': 'CNN_Pose_Score'})
    df_vina = pd.read_csv(os.path.join(INPUT_DIR, "master_vina_affinity.csv")).rename(columns={'Score': 'Vina_Affinity'})

    master_df = df_cnn_aff.merge(df_cnn_pose, on=['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Test_Version'])
    master_df = master_df.merge(df_vina, on=['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Test_Version'])

    master_df = master_df[~master_df['Ligand'].isin(['ATP', 'ADP'])]
    master_df['Vina_Affinity'] = master_df['Vina_Affinity'] * -1
    master_df[['Is_Active', 'NDCG_Weight']] = master_df['Ligand'].apply(assign_ground_truth).apply(pd.Series)

    master_df['Base_Config'] = master_df['Pocket'] + '_' + master_df['Exhaustiveness'].astype(str) + '_' + master_df['Mode'] + '_' + master_df['Hardware']

    # 3. Filter for the Highest Test Version
    master_df['Test_Version'] = master_df['Test_Version'].astype(int)
    idx_max_version = master_df.groupby('Base_Config')['Test_Version'].transform('max') == master_df['Test_Version']
    best_version_df = master_df[idx_max_version].copy()

    results = []

    # 4. Calculate Consensus & Math
    for name, group in best_version_df.groupby(['Base_Config', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware']):
        base_config, pocket, ex, mode, hw = name
        
        group = group.copy()
        group['Norm_CNN_Aff'] = min_max_normalize(group['CNN_Affinity'])
        group['Norm_CNN_Pose'] = min_max_normalize(group['CNN_Pose_Score'])
        group['Norm_Vina'] = min_max_normalize(group['Vina_Affinity'])
        # Calculate Weighted Consensus Score:
        # 25% CNN Affinity, 45% CNN Pose Score, 30% Vina Affinity
        group['Consensus'] = (
            (0.25 * group['Norm_CNN_Aff']) + 
            (0.45 * group['Norm_CNN_Pose']) + 
            (0.30 * group['Norm_Vina'])
        )
        sorted_group = group.sort_values(by='Consensus', ascending=False)
        
        ndcg_score = calculate_ndcg(sorted_group['NDCG_Weight'].tolist())
        rdkit_data = [[flag] for flag in sorted_group['Is_Active'].tolist()]
        bedroc_score = Scoring.CalcBEDROC(rdkit_data, 0, alpha=8.0)

        # CHANGE 1: We only need 'ex' + the number now, since the shape handles the mode.
        label = f"{ex}"

        results.append({
            'Base_Config': base_config,
            'Pocket': pocket,
            'Hardware': hw,
            'Mode': mode,         # ADDED: Keep mode so Seaborn knows what shape to draw
            'Label': label,
            'BEDROC': bedroc_score,
            'NDCG': ndcg_score
        })

    results_df = pd.DataFrame(results)

    # CHANGE 2: Outer Merge
    # We use an 'outer' join so that if a hardware log is missing, the accuracy dot still plots 
    # (it just won't have an X-axis time). We then drop rows where Time_Min is explicitly NaN.
    final_df = results_df.merge(hw_avg, on='Base_Config', how='outer').dropna(subset=['Time_Min', 'BEDROC'])
    
    return final_df

def draw_scatter_plots(df):
    print("Drawing Cost-Benefit Scatterplots...")
    sns.set_theme(style="whitegrid", font="sans-serif")
    
    # CHANGE 3: Only the top row (HEAT) as requested
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    palette = {'GPU': 'forestgreen', 'CPU': 'royalblue'}
    # CHANGE 4: Add the shape dictionary (Circle = rescore, Star = refinement)
    markers = {'rescore': 'o', 'refinement': '*'}
    
    metrics = [('BEDROC', 'Weighted Consensus BEDROC (Alpha=8)'), ('NDCG', 'Consensus NDCG')]

    subset = df[df['Pocket'] == '7B7D_HEAT']
    
    for col_idx, (metric_col, metric_title) in enumerate(metrics):
        ax = axes[col_idx]
        
        if subset.empty:
            ax.set_title("7B7D_HEAT: No Data Found")
            continue
            
        # CHANGE 5: Inject the 'style' and 'markers' arguments into the scatterplot
        sns.scatterplot(
            data=subset, x='Time_Min', y=metric_col, hue='Hardware', 
            style='Mode', markers=markers, palette=palette, 
            s=200, edgecolor='black', ax=ax
        )
        
        # --- THE FIX: Apply Log Scale and Limits BEFORE AdjustText ---
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x:g}"))
        ax.set_ylim(-0.05, 1.1)
        ax.set_title(f"7B7D_HEAT: Compute Time vs {metric_col}", fontsize=14, weight='bold', pad=10)
        ax.set_xlabel("Average Time per Molecule (Minutes)", fontsize=12, weight='bold')
        ax.set_ylabel(metric_title, fontsize=12, weight='bold')

        # --- THE FIX: SMART TEXT LABELING ---
        # Instead of drawing text instantly, we build a list of text objects
        texts = []
        for i, row in subset.iterrows():
            # Fallback Manual Offset: Push GPU up slightly, CPU down slightly
            y_offset = 0.02 if row['Hardware'] == 'GPU' else -0.02
            
            # Create the text but don't finalize its position yet
            txt = ax.text(row['Time_Min'], row[metric_col] + y_offset, row['Label'], 
                          fontsize=10, weight='bold', color='black', ha='center')
            texts.append(txt)

        # Call the physics engine to repel the words from each other
        try:
            from adjustText import adjust_text
            # force_text and force_points dictate how aggressively they repel
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color='gray', lw=0.8, alpha=0.7))
        except ImportError:
            print("\nNotice: 'adjustText' library not found. Basic vertical offsets applied instead.")
            print("To get perfect label separation, run: pip install adjustText\n")
        # ------------------------------------
        
        if col_idx == 1:
            # Tell the legend to display both the colors and the shapes
            ax.legend(title="Settings", loc="lower right", fontsize=11)
        else:
            ax.get_legend().remove()

    plt.tight_layout()
    output_img = os.path.join(OUTPUT_DIR, "Hardware_Mol_time_vs_sorting.png")
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Success! Master visualization saved to {output_img}")

if __name__ == '__main__':
    df_final = process_data()
    if not df_final.empty:
        draw_scatter_plots(df_final)
    else:
        print("Error: DataFrame is empty. Merge failed.")