#!/bin/bash

# Define the list of files to search for
files_to_search=("OUTCAR" "INCAR" "POSCAR" "vasprun.xml")

# Get the name of the current directory and its parent directory
current_dir=$(basename "$PWD")
parent_dir=$(dirname "$PWD")

# Define the destination folder where files will be copied, in the parent directory
destination_dir="${parent_dir}/${current_dir}_collected_files"

# Create the destination directory in the parent directory if it doesn't exist
mkdir -p "$destination_dir"

# Find and copy the desired files while mirroring the directory structure
for file in "${files_to_search[@]}"; do
    find . -name "$file" | while read -r file_path; do
        # Get the directory path of the file
        original_dir=$(dirname "$file_path")
        
        # Create the corresponding directory structure in the destination folder (in parent directory)
        mkdir -p "$destination_dir/$original_dir"
        
        # Copy the file to the mirrored directory in the destination folder
        cp "$file_path" "$destination_dir/$original_dir"
    done
done

# Zip the collected files into a package named after the current directory
zip_file="${parent_dir}/${current_dir}_collected_files.zip"
zip -r "$zip_file" "$destination_dir"

# Clean up the temporary directory (optional)
rm -rf "$destination_dir"

echo "Collected files have been zipped into $zip_file in the parent directory"
