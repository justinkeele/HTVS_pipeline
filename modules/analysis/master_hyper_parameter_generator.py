import pandas as pd
import numpy as np
import glob
import os
import math

print("Calculating advanced virtual screening metrics (Deduplicating CPU/GPU)...")

# ==========================================
# 1. Metric Calculation Functions
# ==========================================

def get_weight(ligand_name):
    """Assigns the ground truth weight based on the filename prefix."""
    if ligand_name.startswith('tophit_high'): return 3
    elif ligand_name.startswith('tophit_low'): return 2
    else: return 0  

def is_hit(ligand_name):
    """Binary classification for Precision/Recall."""
    return 1 if get_weight(ligand_name) > 0 else 0

def calc_precision_at_k(sorted_ligands, k=5):
    top_k = sorted_ligands[:k]
    hits = sum([is_hit(lig) for lig in top_k])
    return hits / k

def calc_recall_at_k(sorted_ligands, total_hits, k=5):
    if total_hits == 0: return 0.0
    top_k = sorted_ligands[:k]
    hits = sum([is_hit(lig) for lig in top_k])
    return hits / total_hits

def calc_dcg(sorted_ligands):
    dcg = 0.0
    for i, lig in enumerate(sorted_ligands):
        rank = i + 1
        weight = get_weight(lig)
        if weight > 0:
            dcg += ( (2**weight - 1) / np.log2(rank + 1) )
    return dcg

def calc_ndcg(sorted_ligands):
    actual_dcg = calc_dcg(sorted_ligands)
    ideal_sorting = sorted(sorted_ligands, key=lambda x: get_weight(x), reverse=True)
    ideal_dcg = calc_dcg(ideal_sorting)
    if ideal_dcg == 0: return 0.0
    return actual_dcg / ideal_dcg

def calc_bedroc(sorted_ligands, alpha=20.0):
    """Calculates BEDROC using the exact RDKit discrete formulation."""
    acts = [is_hit(lig) for lig in sorted_ligands]
    numMol = len(acts)
    numActives = sum(acts)
    if numActives == 0 or numActives == numMol: return 0.0
    
    denom = (1.0 / numMol) * ((1.0 - math.exp(-alpha)) / (math.exp(alpha / numMol) - 1.0))
    sum_exp = 0.0
    
    for i in range(numMol):
        if acts[i]:
            sum_exp += math.exp(-(alpha * (i + 1)) / numMol)
            
    RIE = sum_exp / (numActives * denom)
    ratio = 1.0 * numActives / numMol
    RIEmax = (1.0 - math.exp(-alpha * ratio)) / (ratio * (1.0 - math.exp(-alpha)))
    RIEmin = (1.0 - math.exp(alpha * ratio)) / (ratio * (1.0 - math.exp(alpha)))
    
    if RIEmax != RIEmin:
        return max(0.0, min(1.0, (RIE - RIEmin) / (RIEmax - RIEmin)))
    return 1.0

# ==========================================
# 2. Main Data Extraction Loop
# ==========================================

tsv_files = glob.glob("*/all_ligands_master_summary.tsv")
all_metrics = []

for file in tsv_files:
    raw_run_name = os.path.dirname(file)
    
    # Strip CPU and GPU tags so we can treat them as the exact same run computationally
    clean_run_name = raw_run_name.replace('_GPU', '').replace('_CPU', '')
    
    if '7B7D_HEAT' in clean_run_name: pocket = '7B7D_HEAT'
    elif '7B7D_far' in clean_run_name: pocket = '7B7D_far'
    elif '2IW3' in clean_run_name: pocket = '2IW3'
    else: pocket = 'Other'
    
    df = pd.read_csv(file, sep='\t', names=['Ligand', 'Vina_Affinity', 'CNN_Pose_Score', 'CNN_Affinity'])
    df = df.dropna()
    
    df['Pose_Rank'] = df['CNN_Pose_Score'].rank(ascending=False, method='min')
    df['Affinity_Rank'] = df['CNN_Affinity'].rank(ascending=False, method='min')
    df['Consensus_Rank'] = (df['Pose_Rank'] + df['Affinity_Rank']) / 2.0
    
    cons_list = df.sort_values(by=['Consensus_Rank', 'CNN_Pose_Score'], ascending=[True, False])['Ligand'].tolist()
    pose_list = df.sort_values(by='CNN_Pose_Score', ascending=False)['Ligand'].tolist()
    aff_list  = df.sort_values(by='CNN_Affinity', ascending=False)['Ligand'].tolist()
    vina_list = df.sort_values(by='Vina_Affinity', ascending=True)['Ligand'].tolist()
    
    total_hits = sum([is_hit(lig) for lig in df['Ligand'].tolist()])
    
    for score_type, sorted_list in [('Consensus', cons_list), ('CNN_Pose', pose_list), 
                                    ('CNN_Affinity', aff_list), ('Vina_Affinity', vina_list)]:
        all_metrics.append({
            'Run_Name': clean_run_name,
            'Pocket': pocket,
            'Score_Type': score_type,
            'Precision@5': calc_precision_at_k(sorted_list, k=5),
            'Recall@5': calc_recall_at_k(sorted_list, total_hits, k=5),
            'NDCG': calc_ndcg(sorted_list),
            'BEDROC': calc_bedroc(sorted_list, alpha=20.0)
        })

# ==========================================
# 3. Deduplication and Sorting
# ==========================================

metrics_df = pd.DataFrame(all_metrics)

# Drop duplicates created by reading both _GPU and _CPU folders of the same run
metrics_df = metrics_df.drop_duplicates(subset=['Run_Name', 'Score_Type'])

# Load the human-readable spreadsheet to get the desired Master Sort Order
try:
    human_scores = pd.read_csv('human_readable_scores.csv')
    
    # Clean the human-readable run names just like we did for the metrics
    human_runs_raw = human_scores['Run_Name'].dropna().tolist()
    ordered_clean_runs = []
    
    # Create an ordered, deduplicated list
    for r in human_runs_raw:
        clean = r.replace('_GPU', '').replace('_CPU', '')
        if clean not in ordered_clean_runs:
            ordered_clean_runs.append(clean)
            
    # Apply the categorical sort order
    metrics_df['Run_Name'] = pd.Categorical(metrics_df['Run_Name'], categories=ordered_clean_runs, ordered=True)
    metrics_df = metrics_df.sort_values(by=['Run_Name', 'Score_Type']).dropna(subset=['Run_Name']).reset_index(drop=True)
except FileNotFoundError:
    print("Warning: human_readable_scores.csv not found. Using alphabetical sorting.")
    metrics_df = metrics_df.sort_values(by=['Pocket', 'Run_Name', 'Score_Type']).reset_index(drop=True)

output_filename = "master_hyperparameter.csv"
metrics_df.to_csv(output_filename, index=False)

print(f"Successfully generated metrics and deduplicated runs.")
print(f"Saved to {output_filename}. Ready for graphing.")