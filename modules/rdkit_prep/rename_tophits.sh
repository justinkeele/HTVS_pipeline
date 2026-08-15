#!/bin/bash


# Define the folder and the name map file
TARGET_DIR="72pH_meeko_ready"


echo "Renaming tophits..."

# Loop through every line in the text file new_name
while IFS= read -r new_name; do
    
    # Skip any empty lines just in case
    [ -z "$new_name" ] && continue

    # Use grep to extract just the "Z123..." part from the line
    Z_id=$(echo "$new_name" | grep -o "Z[0-9]*")

    # If the Z-ID exists AND the corresponding .sdf file exists, rename it
    if [ -n "$Z_id" ] && [ -f "$TARGET_DIR/${Z_id}.sdf" ]; then
        mv "$TARGET_DIR/${Z_id}.sdf" "$TARGET_DIR/${new_name}.sdf"
        echo "Renamed: ${Z_id}.sdf -> ${new_name}.sdf"
    fi

done < molecule_names.txt

echo ""
echo "Maybe it worked!"
