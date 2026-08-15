#!/bin/bash

# Loop through the scored mol2 files
for file in *_docked_scored.mol2; do

    # Grep the score line, and use awk to grab just the 3rd column (the number)
    score=$(grep "Grid_Score:" "$file" | awk '{print $3}')

    # Strip the messy suffix off the filename so the table is readable
    clean_name=$(basename "$file" _docked_scored.mol2)

    if [ ! -z "$score" ]; then
        # Print the score first, then the name, separated by a tab
        echo -e "$score\t$clean_name"
    fi

# Pipe the entire loop output into sort. 
# -n sorts numerically (most negative at the top)
done | sort -n > sorted_docking_scores.txt