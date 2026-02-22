#!/bin/bash

# Directory containing the files
directory=$1

# Loop through all matching files in the directory
for file in "$directory"/*rhi1.vert.mmclx.gz; do
    # Check if the file exists to avoid issues when no files match
    if [[ -f "$file" ]]; then
	# Construct the new filename
	new_file="${file/.rhi1.vert.mmclx.gz/.rhi1.mmclx.gz}"
	# Rename the file
	mv "$file" "$new_file"
	echo "Renamed: $file -> $new_file"
    fi
done

echo "Renaming complete!"
