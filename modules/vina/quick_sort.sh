#!/bin/bash

for file in *_results.log; do

    score=$(sed -n '14p' "$file" | awk '{print $5}')

    if [ ! -z "$score" ]; then
        echo -e "$file\t$score"
    fi

done| sort -k2,2nr > test_sorted_docking_scores.txt
