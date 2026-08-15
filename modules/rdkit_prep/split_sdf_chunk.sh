#!/bin/bash

awk '
{
    # 1. Add every line we read into a temporary holding tank
    buffer = buffer $0 "\n"
}
/>  <CatalogID>/ {
    # 2. When we see the CatalogID tag, grab the very next line
    getline id
    buffer = buffer id "\n"
    
    # 3. Strip out any invisible Windows carriage returns or weird characters
    name = id
    gsub(/[^a-zA-Z0-9_]/, "", name)
    
    # 4. Set the name of the file we are about to create
    filename = name ".sdf"
}
/^\$\$\$\$/ {
    # 5. When we hit the $$$$ delimiter, dump the holding tank into the new file
    if (filename != ".sdf" && filename != "") {
        printf "%s", buffer > filename
    }
    
    # 6. Empty the tank and clear the filename for the next molecule
    buffer = ""
    filename = ""
}' chunk_1.sdf

#Add the correct filename you wish to split above.
