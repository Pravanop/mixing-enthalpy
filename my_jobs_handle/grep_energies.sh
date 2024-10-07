#!/bin/bash

# Get the name of the current directory
current_dir=$(basename "$PWD")

# Initialize lists to store results
wavefunctions_found=()
wavefunctions_not_found=()

# Output JSON file named after the current directory
json_output_file="${current_dir}_energy.json"
echo "{" > "$json_output_file"  # Initialize the JSON object

# Loop over all directories inside the 'relaxed' folder
for dir in relaxed/*/; do
    # Define the path to the 'stdout' and 'OUTCAR' files
    static_stdout="${dir}static/stdout"
    outcar_file="${dir}OUTCAR"
    
    # Extract the base directory name (e.g., "V" from "relaxed/V/")
    base_dir=$(basename "$dir")

    # Check if the stdout file exists in the static subdirectory
    if [ -f "$static_stdout" ]; then
        # Search for the pattern "writing wavefunctions" in the stdout file
        if grep -q "writing wavefunctions" "$static_stdout"; then
            # If found, add the directory to the 'wavefunctions_found' list
            wavefunctions_found+=("${base_dir}")
            
            # Now check the OUTCAR file for the "free  energy" pattern
            if [ -f "$outcar_file" ]; then
                # Extract the line with "free  energy" and store the energy value
                energy_line=$(grep "free  energy" "$outcar_file")
                if [[ $energy_line =~ TOTEN[[:space:]]+=[[:space:]]+([-0-9.]+) ]]; then
                    energy_value="${BASH_REMATCH[1]}"

                    # Now extract the number of ions (NIONS) from OUTCAR
                    nions_line=$(grep "NIONS" "$outcar_file")
                    if [[ $nions_line =~ NIONS[[:space:]]+=[[:space:]]+([0-9]+) ]]; then
                        nions="${BASH_REMATCH[1]}"

                        # Calculate energy per ion
                        energy_per_ion=$(echo "scale=8; $energy_value / $nions" | bc)

                        # Add the directory name and energy per ion to the JSON file
                        echo "\"${base_dir}\": \"${energy_per_ion}\"," >> "$json_output_file"
                    else
                        echo "Error: Could not extract NIONS from ${outcar_file}"
                    fi
                else
                    echo "Error: Could not extract energy value from ${outcar_file}"
                fi
            else
                echo "Warning: OUTCAR file not found in ${dir}"
            fi
        else
            # If "writing wavefunctions" is not found, add to the not found list
            wavefunctions_not_found+=("${base_dir}")
        fi
    else
        echo "Warning: stdout file not found in ${dir}static"
    fi
done

# Close the JSON object (remove the trailing comma)
sed -i '$ s/,$//' "$json_output_file"
echo "}" >> "$json_output_file"

# Display the results
echo "Directories where 'writing wavefunctions' was found in static/stdout:"
printf "%s\n" "${wavefunctions_found[@]}"

echo -e "\nDirectories where 'writing wavefunctions' was NOT found in static/stdout:"
printf "%s\n" "${wavefunctions_not_found[@]}"

echo -e "\nEnergy per ion values have been written to $json_output_file."
