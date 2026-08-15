#!/bin/bash

# Temporary workspace to hold all 3 scores on one line
tmp_file="master_scores.tmp"
> "$tmp_file"

for file in *_results.log; do
    # Extract Mode 1: $2 = Vina (kcal/mol), $4 = CNNpose (0-1), $5 = CNNaffinity (pK)
    scores=$(awk '/^ *1 / {print $2, $4, $5; exit}' "$file")

    if [ -n "$scores" ]; then
        # Strip "_results.log" so the table only shows clean molecule names
        clean_name=$(basename "$file" _results.log)
        echo -e "${clean_name}\t${scores}" >> "$tmp_file"
    fi
done

# 1. Vina / Vinardo score (Column 2 -> ascending -n so most negative is top)
sort -k2,2n "$tmp_file" | awk '{print $1 "\t" $2}' > sorted_vina_affinity.txt

# 2. CNN Pose score (Column 3 -> descending -nr so closest to 1.0 is top)
sort -k3,3nr "$tmp_file" | awk '{print $1 "\t" $3}' > sorted_cnn_pose_score.txt

# 3. CNN Affinity (Column 4 -> descending -nr so highest pK is top)
sort -k4,4nr "$tmp_file" | awk '{print $1 "\t" $4}' > sorted_cnn_affinity.txt

#=====================================================================
# Consensus Filtering & CNN_VS Ranking
#=====================================================================
# Parameters Justification:
#   - Vina ($2) <= -6.0 kcal/mol : Thermodynamic hard gate against small decoys
#   - CNNpose ($3) >= 0.50       : Koes lab standard threshold for <2A RMSD confidence
#   - CNN_VS ($3 * $4)           : Weighted product of pose confidence * predicted pK

awk '
{
    vina = $2
    cnn_pose = $3
    cnn_aff = $4
    
    # Apply Thermodynamic Gate AND Neural Network Confidence Gate
    if (vina <= -6.0 && cnn_pose >= 0.50) {
        cnn_vs = cnn_pose * cnn_aff
        printf "%s\t%.4f\t(Vina: %.2f | Pose: %.4f | pK: %.3f)\n", $1, cnn_vs, vina, cnn_pose, cnn_aff
    }
}' "$tmp_file" | sort -k2,2nr > sorted_consensus_score.txt

# Clean up the temp file
rm -f "$tmp_file"

echo "Done! Created 3 sorted files:"
echo "  1. sorted_vina_affinity.txt   (Empirical Vina energy - kcal/mol)"
echo "  2. sorted_cnn_pose_score.txt  (CNN binding mode probability - 0.0 to 1.0)"
echo "  3. sorted_cnn_affinity.txt    (CNN predicted pK affinity)"
