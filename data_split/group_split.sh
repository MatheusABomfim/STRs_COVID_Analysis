#!/bin/bash

# Directory where the files are located (defined directly)
input_dir=~/projects/strling/STRling/cbgm_2025/cbgm_output

# Output directory (defined directly)
output_dir=/home/matheusbomfim/projects/strs_paper/samples

# Verify if input directory exists
if [[ ! -d "$input_dir" ]]; then
  echo "Input directory '$input_dir' does not exist."
  exit 1
fi

# Create output directories if they do not exist
mkdir -p "$output_dir/case" "$output_dir/control" "$output_dir/gp_global"

# Group CSV file
group_csv="grupos.csv"

# Check if the grupos.csv file exists
if [[ ! -f "$group_csv" ]]; then
  echo "File $group_csv not found!"
  exit 1
fi

# Arrays to store samples by group
declare -A groups  # Associative array: groups["name"]="group"

# Read CSV filling the array
while IFS=, read -r sample group; do
  # Remove spaces and line breaks
  sample=$(echo "$sample" | tr -d '[:space:]')
  group=$(echo "$group" | tr -d '[:space:]')
  
  # Ignore header and empty lines
  [[ "$sample" == "sample" || -z "$sample" ]] && continue
  
  # Store in array (keeping .bam if it exists)
  groups["$sample"]="$group"
done < "$group_csv"

echo "Total samples loaded: ${#groups[@]}"
echo "  - Case: $(grep -c ",case" "$group_csv")"
echo "  - Control: $(grep -c ",control" "$group_csv")"

# Counters  
count_case=0
count_control=0
count_none=0

# Process genotype files
for file in "$input_dir"/*-genotype.txt; do
  filename=$(basename "$file")
  
  # Extract sample_id keeping consistency with CSV
  sample_id=$(echo "$filename" | sed 's/-genotype\.txt$//')
  
  # Check if the sample is in the CSV
  if [[ -n "${groups[$sample_id]}" ]]; then
    group="${groups[$sample_id]}"
    
    # Copy files according to the group
    if [[ "$group" == "case" ]]; then
      cp -v "$file" "$output_dir/case/"
      ((count_case++))
    elif [[ "$group" == "control" ]]; then
      cp -v "$file" "$output_dir/control/"
      ((count_control++))
    fi
    
    # Copy to gp_global in both cases
    cp -v "$file" "$output_dir/gp_global/"
    
    # Copy corresponding unplaced file
    unplaced_file="$input_dir/${sample_id}-unplaced.txt"
    if [[ -f "$unplaced_file" ]]; then
      cp -v "$unplaced_file" "$output_dir/gp_global/"
    else
      echo "Warning: Unplaced file not found for $sample_id"
    fi
    
  else
    echo "File $filename does not belong to any group"
    ((count_none++))
  fi
done

# Summary
echo
echo "===== Summary ====="
echo "Samples copied to 'case': $count_case"
echo "Samples copied to 'control': $count_control"
echo "Samples copied to gp_global: $((count_case + count_control))"
echo "Samples NOT included: $count_none"

# Final check
echo
echo "Checking copied files:"
echo "Case: $(ls -1 "$output_dir/case/" | wc -l)"
echo "Control: $(ls -1 "$output_dir/control/" | wc -l)"
echo "Global: $(ls -1 "$output_dir/gp_global/" | wc -l)"