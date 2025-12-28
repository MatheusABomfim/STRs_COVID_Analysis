#!/bin/bash

# ============================================================================
# STR DATA CONSOLIDATION SCRIPT - Processes outliers_gp_global*.STRs.tsv files
# ============================================================================
# This script consolidates multiple STR outliers files into a single TSV file
# and generates a statistical summary per sample.
# ============================================================================


# Configurações
dir_vcfs="../samples/gp_global"
output_file="${dir_vcfs}/merged/outliers_consolidated.tsv"
sample_summary="${dir_vcfs}/merged/outliers_summary.tsv"
log_file="${dir_vcfs}/merged/outliers_merge_log.txt"

# Create directory 
mkdir -p "$(dirname "$output_file")"

# Iniciate log file
{
echo "======================================"
echo "STR DATA CONSOLIDATION (OUTLIERS FILES)"
echo "Started: $(date)"
echo "======================================"
echo "Input directory: $dir_vcfs"
echo "Output file: $output_file"
echo "Sample summary: $sample_summary"
echo "Log file: $log_file"
echo "--------------------------------------"
} | tee "$log_file"

# Verify with directoy exists
if [ ! -d "$dir_vcfs" ]; then
    echo "ERROR: Directory $dir_vcfs does not exist." | tee -a "$log_file"
    exit 1
fi

# search outliers_gp_global* files
shopt -s nullglob
TSV_FILES=("$dir_vcfs"/outliers_gp_global*.STRs.tsv)

# Verify if the files was identified
if [ ${#TSV_FILES[@]} -eq 0 ]; then
    echo "ERROR: No outliers_gp_global*.STRs.tsv files found in $dir_vcfs" | tee -a "$log_file"
    echo "Available files:" | tee -a "$log_file"
    ls -la "$dir_vcfs"/*.tsv 2>/dev/null | head -20 | tee -a "$log_file"
    exit 1
fi

echo "Found ${#TSV_FILES[@]} outliers files" | tee -a "$log_file"
echo "--------------------------------------" | tee -a "$log_file"

# Create headers of consolidate file
echo -e "sample\tchrom\tstart\tend\trepeat_unit\tallele1_est\tallele2_est\tspanning_reads\tspanning_pairs\tleft_clips\tright_clips\tunplaced_pairs\tsum_str_counts\tdepth\toutlier\tp\tp_adj" > "$output_file"
echo -e "sample\tvariant_count\tmean_depth\tmean_allele_diff\thomozygous_count\theterozygous_count" > "$sample_summary"

# Counters
counter=0
total_variants=0
total_samples=${#TSV_FILES[@]}
failed_samples=0

# Process every file
for TSV_FILE in "${TSV_FILES[@]}"; do
    counter=$((counter + 1))
    
    # Extract nome of the samples
    BASE_NAME=$(basename "$TSV_FILE" .STRs.tsv)
    SAMPLE_ID=${BASE_NAME#outliers_gp_global}
    
    # Exclude the extension (.bam)
    SAMPLE_ID=${SAMPLE_ID%.bam}
    
    echo -n "[$counter/$total_samples] Processing $SAMPLE_ID from $(basename "$TSV_FILE")... " | tee -a "$log_file"
    
    # Verify if the file is empty
    if [ ! -s "$TSV_FILE" ]; then
        echo "SKIPPED - File is empty" | tee -a "$log_file"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Verify the rows along the file
    line_count=$(wc -l < "$TSV_FILE")
    if [ "$line_count" -le 1 ]; then
        echo "SKIPPED - File has only header or no data" | tee -a "$log_file"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Process the files with awk
    awk_result=$(awk -v sample="$SAMPLE_ID" '
    BEGIN {
        FS = "\t"
        OFS = "\t"
        total = 0
        depth_sum = 0
        diff_sum = 0
        hom = 0
        het = 0
        errors = 0
        
        # Map fields
        # chrom   left    right   locus   sample  group   repeatunit      allele1_est     allele2_est     spanning_reads  spanning_pairs  left_clips      right_clips     unplaced_pairs  sum_str_counts  sum_str_log   depth   outlier p       p_adj
    }
    
    # Jump lines 
    NR == 1 && /^chrom.*/ { next }
    /^#/ { next }
    
    # Process data line by line
    NF >= 17 {
        # Basic fields
        chrom = $1
        start = $2
        end = $3
        repeat_unit = $7
        outlier = $18
        p_value = $19
        p_adj = $20
        
        # Manage the alelles
        allele1_str = $8
        allele2_str = $9
        
        if (allele1_str == "NaN" || allele1_str == "nan" || allele1_str == "" || allele1_str == ".") {
            allele1 = 0
        } else {
            allele1 = allele1_str + 0
        }
        
        if (allele2_str == "NaN" || allele2_str == "nan" || allele2_str == "" || allele2_str == ".") {
            allele2 = 0
        } else {
            allele2 = allele2_str + 0
        }
        
        # Numeric values
        spanning_reads = $10 + 0
        spanning_pairs = $11 + 0
        left_clips = $12 + 0
        right_clips = $13 + 0
        unplaced_pairs = $14 + 0
        sum_str_counts = $15 + 0
        
        # Depth 
        depth_field = 17
        if (NF == 19) { 
            depth_field = 16
            outlier = $17
            p_value = $18
            p_adj = $19
        }
        depth = $(depth_field) + 0
        
        # Validate data
        if (chrom == "" || start == "" || end == "") {
            errors++
            next
        }
        
        # Statistics
        total++
        depth_sum += depth
        
        # Calculate the difference between alelles
        diff = (allele2 > allele1) ? allele2 - allele1 : allele1 - allele2
        diff_sum += diff
        
        # Count homozygous vs hetezygous
        if (int(allele1 + 0.5) == int(allele2 + 0.5)) {
            hom++
        } else {
            het++
        }
        
        # Print the consolidate header
        print sample, chrom, start, end, repeat_unit, allele1, allele2, \
              spanning_reads, spanning_pairs, left_clips, right_clips, \
              unplaced_pairs, sum_str_counts, depth, outlier, p_value, p_adj
    }
    
    END {
        if (total > 0) {
            mean_depth = depth_sum / total
            mean_diff = diff_sum / total
            printf "%s\t%d\t%.2f\t%.2f\t%d\t%d\n", sample, total, mean_depth, mean_diff, hom, het
        } else if (errors > 0) {
            printf "%s\t0\t0.00\t0.00\t0\t0\n", sample
        }
    }
    ' "$TSV_FILE")
    
    # Extract statistics
    stats_line=$(echo "$awk_result" | tail -1)
    variant_data=$(echo "$awk_result" | head -n -1 2>/dev/null || echo "")
    
    # Add data to final file
    if [ -n "$variant_data" ]; then
        echo "$variant_data" >> "$output_file"
    fi
    
    # Add statitics to summary
    if [ -n "$stats_line" ] && echo "$stats_line" | grep -q $'\t'; then
        echo "$stats_line" >> "$sample_summary"
        sample_variants=$(echo "$stats_line" | cut -f2)
        total_variants=$((total_variants + sample_variants))
    fi
    
    # Count variants per file
    variant_count_num=$(echo "$variant_data" | wc -l 2>/dev/null || echo "0")
    echo "OK - $variant_count_num variants" | tee -a "$log_file"
done

# Finish the log
{
echo "======================================"
echo "CONSOLIDATION COMPLETE"
echo "======================================"
echo "SUMMARY STATISTICS:"
echo "  Total samples found: $total_samples"
echo "  Successfully processed: $((total_samples - failed_samples))"
echo "  Failed/skipped: $failed_samples"
echo "  Total variants consolidated: $total_variants"
echo ""

if [ $((total_samples - failed_samples)) -gt 0 ]; then
    avg_variants=$((total_variants / (total_samples - failed_samples)))
    echo "  Average variants per sample: $avg_variants"
fi

echo ""
echo "OUTPUT FILES:"
if [ -f "$output_file" ]; then
    output_size=$(du -h "$output_file" | cut -f1)
    output_lines=$(wc -l < "$output_file")
    echo "  Consolidated TSV: $output_file"
    echo "    Size: $output_size"
    echo "    Lines: $((output_lines - 1)) variants + 1 header"
else
    echo "  ERROR: Consolidated TSV was not created"
fi

if [ -f "$sample_summary" ]; then
    summary_size=$(du -h "$sample_summary" | cut -f1)
    summary_lines=$(wc -l < "$sample_summary")
    echo "  Sample summary: $sample_summary"
    echo "    Size: $summary_size"
    echo "    Lines: $((summary_lines - 1)) samples + 1 header"
else
    echo "  ERROR: Sample summary was not created"
fi

echo ""
echo "First 3 lines of consolidated output:"
head -n 4 "$output_file" 2>/dev/null || echo "  (File not available)"

echo ""
echo "First 3 samples in summary:"
head -n 4 "$sample_summary" 2>/dev/null || echo "  (File not available)"

echo ""
echo "Finished: $(date)"
echo "======================================"
} | tee -a "$log_file"

# Copy log to stdout
echo ""
echo "Log saved to: $log_file"
cat "$log_file" | tail -20

exit 0
