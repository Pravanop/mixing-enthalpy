#!/bin/bash

# Initialize lists to store directories
found_list=()
relaxed_list=()

# Create a 'relaxed' directory if it doesn't exist
if [ ! -d "relaxed" ]; then
    mkdir relaxed
fi

# Loop over all stdout files in subdirectories
for dir in */; do
    stdout_file="${dir}stdout"
    
    # Check if stdout file exists
    if [ ! -f "$stdout_file" ]; then
        echo "Warning: stdout file not found in ${dir}"
        continue
    fi

    # Check if the pattern "reached required accuracy" is found first
    if grep -q "reached required accuracy" "$stdout_file"; then
        # If the accuracy pattern is found, consider the job done, move it to relaxed
        relaxed_list+=("${dir}")
        mv "$dir" relaxed/
        
        # Navigate into the 'static' subdirectory, copy CONTCAR to POSCAR, and run sbatch
        if [ -d "relaxed/${dir}static" ]; then
            (
                cd "relaxed/${dir}static" || exit  # Navigate into the 'static' subdirectory
                cp ../CONTCAR ./POSCAR  # Copy CONTCAR to POSCAR
                sbatch runjob  # Submit the job with sbatch (commented out as per your modification)
            )
        else
            echo "Warning: static directory not found in relaxed/${dir}"
        fi
    # If "reached required accuracy" is not found, check for the "rerun" pattern
    elif grep -q "rerun" "$stdout_file"; then
        # If rerun is found, re-execute the commands
        found_list+=("${dir}")
        
        # Navigate to the directory and run the required commands
        (
            cd "$dir" || exit  # Navigate into the directory
            cp CONTCAR POSCAR  # Copy CONTCAR to POSCAR
            sbatch runjob  # Submit the job with sbatch (commented out as per your modification)
        )
    fi
done

# Display the lists
echo "Directories where pattern 'rerun' was found and commands executed:"
printf "%s\n" "${found_list[@]}"

echo -e "\nDirectories where pattern 'reached required accuracy' was found, moved to 'relaxed', and sbatch runjob executed in 'static':"
printf "%s\n" "${relaxed_list[@]}"

