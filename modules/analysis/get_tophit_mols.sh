#!/bin/bash


# The mega sort
mega_file=mega_sorted_docking_scores.txt
top_hits=top_hits.txt

find . -name "sorted_docking_scores.txt" -exec cat {} + | sort -k2,2n > "$mega_file"
echo "sort completed"
wc -l "$mega_file" | awk '{print "Total lines in file: " $1}'

head -n 200 "$mega_file" > top_hits.txt

# Get our pdbqt and log files to look at
echo "getting pdbqt and log files for top hits"

mkdir -p tophit_pdbqts
mkdir -p tophit_logs

awk '{print $1}' "$top_hits" | while read -r zinc_id; do
    cp batch_*/"${zinc_id}.pdbqt" tophit_pdbqts/ 2>/dev/null
    cp batch_*/"${zinc_id}_results.log" tophit_logs/ 2>/dev/null
done

pdbqt_count=$(ls -1 tophit_pdbqts/*.pdbqt 2>/dev/null | wc -l)
log_count=$(ls -1 tophit_logs/*_results.log 2>/dev/null | wc -l)

echo "pdbqt files copied: $pdbqt_count"
echo "log files copied: $log_count"