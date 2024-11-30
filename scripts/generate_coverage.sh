#!/bin/bash

# Paths for config files
COVERAGE_YAML="config/coverage.yaml"
CONFIG_YAML="config/config.yaml"

# Remove any existing coverage.yaml file and start fresh
rm -f "$COVERAGE_YAML"
touch "$COVERAGE_YAML"

# Initialize the genome_coverage section in coverage.yaml
echo -e "\ngenome_coverage: {" >> "$COVERAGE_YAML"

# Initialize empty arrays for shallow and deep samples
shallow_samples=()
deep_samples=()

# Loop through each sample directory and add to genome_coverage
for sample_dir in "${PWD}/13_flyerough"/*; do
    sample=$(basename "$sample_dir")
    depthcheck_file="$sample_dir/flye01/depthcheckflag.txt"
    
    # Check if the depthcheckflag.txt file exists and read its content
    if [[ -f "$depthcheck_file" ]]; then
        coverage=$(cat "$depthcheck_file")
        echo "  \"$sample\": \"$coverage\"," >> "$COVERAGE_YAML"

        # Check if coverage is less than 70 and add to shallow_samples if true
        if (( coverage < 70 )); then
            shallow_samples+=("$sample")
        else
            deep_samples+=("$sample")
        fi
    fi
done

# Remove the last comma and add closing bracket for genome_coverage
sed -i '$ s/,$//' "$COVERAGE_YAML"
echo "}" >> "$COVERAGE_YAML"

# Append shallow_sample_names section to coverage.yaml
echo -e "\nshallow_sample_names: [" >> "$COVERAGE_YAML"
for sample in "${shallow_samples[@]}"; do
    echo "  \"$sample\"," >> "$COVERAGE_YAML"
done
sed -i '$ s/,$//' "$COVERAGE_YAML"  # Remove the last comma
echo "]" >> "$COVERAGE_YAML"        # Close shallow_sample_names list

# Append deep_sample_names section to coverage.yaml
echo -e "\ndeep_sample_names: [" >> "$COVERAGE_YAML"
for sample in "${deep_samples[@]}"; do
    echo "  \"$sample\"," >> "$COVERAGE_YAML"
done
sed -i '$ s/,$//' "$COVERAGE_YAML"  # Remove the last comma
echo "]" >> "$COVERAGE_YAML"        # Close deep_sample_names list

# Copy the final content of coverage.yaml to config.yaml
cat "$COVERAGE_YAML" >> "$CONFIG_YAML"

# Create a flag file to indicate the script has run successfully
touch "${PWD}/13_flyerough/updateconfigflag.txt"
