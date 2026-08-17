import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Configuration & Paths
INPUT_DIR = "../../../docs/graphs_with_good_paths"
OUTPUT_DIR = "../../../docs/graphs_with_good_paths"

X_AXIS_ORDER = [
    'rescore_ex8', 'rescore_ex16', 'rescore_ex32', 'rescore_ex64', 'rescore_ex128', 
    'refinement_ex8', 'refinement_ex16', 'refinement_ex32'
]

# We define the categories for our color mapping
def get_ligand_category(ligand_name):
    if "tophit_high" in ligand_name:
        return "tophit_high"
    elif "tophit_low" in ligand_name:
        return "tophit_low"
    elif "no_inhibiton" in ligand_name:
        return "no_inhibition"
    else:
        return "z_series"

def min_max_normalize(series):
    '''Squashes scores strictly between 0.0 and 1.0.'''
    if series.max() == series.min(): 
        return series * 0.0 + 0.5
    return (series - series.min()) / (series.max() - series.min())

def process_data():
    print("Loading master CSVs...")
    df_cnn_aff = pd.read_csv(os.path.join(INPUT_DIR, "master_cnn_affinity.csv")).rename(columns={'Score': 'CNN_Affinity'})
    df_cnn_pose = pd.read_csv(os.path.join(INPUT_DIR, "master_cnn_pose_score.csv")).rename(columns={'Score': 'CNN_Pose_Score'})
    df_vina = pd.read_csv(os.path.join(INPUT_DIR, "master_vina_affinity.csv")).rename(columns={'Score': 'Vina_Affinity'})

    master_df = df_cnn_aff.merge(df_cnn_pose, on=['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Test_Version'])
    master_df = master_df.merge(df_vina, on=['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Test_Version'])

    # Drop ATP and ADP
    master_df = master_df[~master_df['Ligand'].isin(['ATP', 'ADP'])]

    # Standardize Vina Affinity direction so higher is better
    master_df['Vina_Affinity'] = master_df['Vina_Affinity'] * -1

    # Map the color categories
    master_df['Category'] = master_df['Ligand'].apply(get_ligand_category)
    master_df['Test_Version_int'] = master_df['Test_Version'].astype(int)
    master_df['X_Label'] = master_df['Mode'] + '_ex' + master_df['Exhaustiveness'].astype(str)

    # --- THE FIX ---
    # Removed the extra '_ex' string because master_df['Exhaustiveness'] already contains it
    master_df['X_Label'] = master_df['Mode'] + '_' + master_df['Exhaustiveness'].astype(str)

    # 2. FILTERING GATE: Remove Test_Version 1 for HEAT and Far pockets
    invalid_mask = (master_df['Pocket'].isin(['7B7D_HEAT', '7B7D_far'])) & (master_df['Test_Version_int'] == 1)
    valid_df = master_df[~invalid_mask].copy()

    # 3. GET HIGHEST VERSION: Group by hardware setup and keep max Test_Version
    valid_df['Config_Key'] = valid_df['Pocket'] + '_' + valid_df['X_Label'] + '_' + valid_df['Hardware']
    idx_max_version = valid_df.groupby('Config_Key')['Test_Version_int'].transform('max') == valid_df['Test_Version_int']
    best_runs_df = valid_df[idx_max_version].copy()

    processed_chunks = []

    # 4. APPLY MATH: Normalize and calculate Consensus exactly like the heatmap
    for name, group in best_runs_df.groupby('Config_Key'):
        if group['X_Label'].iloc[0] not in X_AXIS_ORDER:
            continue
            
        group = group.copy()
        group['Norm_CNN_Aff'] = min_max_normalize(group['CNN_Affinity'])
        group['Norm_CNN_Pose'] = min_max_normalize(group['CNN_Pose_Score'])
        group['Norm_Vina'] = min_max_normalize(group['Vina_Affinity'])
        
        # Weighted Consensus: 25% Aff, 45% Pose, 30% Vina
        group['Consensus'] = (
            (0.25 * group['Norm_CNN_Aff']) + 
            (0.45 * group['Norm_CNN_Pose']) + 
            (0.30 * group['Norm_Vina'])
        )
        
        processed_chunks.append(group)

    return pd.concat(processed_chunks, ignore_index=True)

def draw_stripplots(df):
    print("Generating Strip-Plot Grid...")
    
    sns.set_theme(style="whitegrid")
    
    metrics = ['Consensus', 'CNN_Affinity', 'CNN_Pose_Score', 'Vina_Affinity']
    pockets = ['7B7D_HEAT', '7B7D_far', '2IW3']
    
    # Map your exact requested colors to the categories
    # Used 'gold' for yellow as it displays much better on white backgrounds
    palette = {
        'tophit_high': 'forestgreen', 
        'tophit_low': 'royalblue', 
        'no_inhibition': 'gold', 
        'z_series': 'crimson'
    }

    # sharey='row' gives each metric its own scale, but links the 3 pockets together
    fig, axes = plt.subplots(len(metrics), len(pockets), figsize=(18, 24), sharex=True, sharey='row')

    for r_idx, metric in enumerate(metrics):
        for c_idx, pocket in enumerate(pockets):
            ax = axes[r_idx, c_idx]
            
            subset = df[df['Pocket'] == pocket]
            
            # --- THE FIX: REMOVE VINA REFINEMENT CLASHES ---
            # If we are graphing the Vina_Affinity row, drop all refinement runs
            # so the extreme steric clash penalties don't blow out the Y-axis scale.
            if metric == 'Vina_Affinity':
                subset = subset[subset['Mode'] == 'rescore']
            
            if subset.empty:
                ax.set_title(f"{pocket} - No Data", fontsize=12)
                continue
                
            # Draw the Stripplot
            # jitter=True prevents the dots from stacking perfectly on top of each other
            sns.stripplot(
                data=subset, x='X_Label', y=metric, hue='Category', 
                palette=palette, order=X_AXIS_ORDER, jitter=0.25, 
                alpha=0.85, size=7, edgecolor='black', linewidth=0.5, ax=ax
            )
            
            # The Vertical Splitter Line Trick (Between index 4 and 5)
            ax.axvline(x=4.5, color='gray', linestyle='--', linewidth=2, alpha=0.7)

            # Formatting
            if r_idx == 0:
                ax.set_title(pocket, fontsize=16, weight='bold', pad=15)
            
            if c_idx == 0:
                # Clean up the Y-axis names for the graph
                clean_metric = metric.replace('_', ' ')
                ax.set_ylabel(clean_metric, fontsize=14, weight='bold')
            else:
                ax.set_ylabel("")
            
            ax.set_xlabel("")
            ax.tick_params(axis='x', rotation=45)
            
            # Handle the legend: We only want one master legend for the whole figure
            if r_idx == 0 and c_idx == 2:
                ax.legend(title="Ligand Category", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
            else:
                if ax.get_legend() is not None:
                    ax.get_legend().remove()

            # --- THE FIX: DYNAMIC MINOR GRIDLINES ---
            import matplotlib.ticker as ticker
            
            # Turn on the minor ticks for the Y-axis
            ax.minorticks_on()
            # We don't want minor ticks on the X-axis because it's categorical text
            ax.xaxis.set_tick_params(which='minor', bottom=False)

            if metric == 'Vina_Affinity':
                # Vina scores range from ~0 to ~15.
                # Put a major line every 2 points, and a minor dashed line every 1 point.
                ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
                ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
            else:
                # CNN and Consensus scores range from 0.0 to 1.0.
                # Put a major line every 0.2, and a minor dashed line every 0.1.
                ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
                ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            
            # Draw the physical lines onto the graph
            ax.grid(True, which='major', axis='y', color='gray', linestyle='-', alpha=0.4)
            ax.grid(True, which='minor', axis='y', color='lightgray', linestyle='--', alpha=0.5)
            # ----------------------------------------
    
    plt.tight_layout()
    output_img = os.path.join(OUTPUT_DIR, "Raw_Scores_Stripplots.png")
    plt.savefig(output_img, dpi=400, bbox_inches='tight')
    print(f"Success! Master Stripplot saved to {output_img}")

if __name__ == '__main__':
    df_results = process_data()
    if not df_results.empty:
        draw_stripplots(df_results)
    else:
        print("Error: No data successfully processed.")