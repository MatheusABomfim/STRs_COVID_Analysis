#!/bin/bash

dir_vcfs="../samples/gp_global"
output_file="${dir_vcfs}/merged/str_consolidated.tsv"  
sample_summary="${dir_vcfs}/merged/sample_summary.tsv"

mkdir -p "$(dirname "$output_file")"
if [ ! -d "$dir_vcfs" ]; then
    echo "ERROR: Directory $dir_vcfs does not exist."
    exit 1
fi

shopt -s nullglob
TSV_FILES=("$dir_vcfs"/*-genotype.txt)
if [ ${#TSV_FILES[@]} -eq 0 ]; then
    echo "No *-genotype.txt files"
    exit 1
fi

echo "======================================"
echo "STR DATA CONSOLIDATION"
echo "======================================"
echo "Input directory: $dir_vcfs"
echo "Found ${#TSV_FILES[@]} genotype files"
echo "Output file: $output_file"
echo "--------------------------------------"

echo "Creating consolidated TSV header..."

# Generate unified file
echo -e "sample\tchrom\tstart\tend\trepeat_unit\tallele1_est\tallele2_est\tanchored_reads\tspanning_reads\tspanning_pairs\texpected_spanning_pairs\tspanning_pairs_pctl\tleft_clips\tright_clips\tunplaced_pairs\tdepth" > "$output_file"

echo -e "sample\tvariant_count\tmean_depth\tmean_allele_diff\thomozygous_count\theterozygous_count" > "$sample_summary"

counter=0
total_variants=0

for TSV_FILE in "${TSV_FILES[@]}"; do
    counter=$((counter + 1))
    BASE_NAME=$(basename "$TSV_FILE" .txt)
    SAMPLE_ID=$(echo "$BASE_NAME" | sed 's/\.bam-genotype$//')
    
    echo -n "[$counter/${#TSV_FILES[@]}] Processing $SAMPLE_ID... "
    
    # Process with awk to handle NaN and calculate stats
    awk_result=$(awk -v sample="$SAMPLE_ID" '
    BEGIN {
	OFS="\t"
        total=0
        depth_sum=0
        diff_sum=0
        hom=0
        het=0
    }
    NR>1 {
        # Substituir NaN por 0
        allele1 = $5
        allele2 = $6
        if (allele1 == "NaN" || allele1 == "nan") allele1 = 0
        if (allele2 == "NaN" || allele2 == "nan") allele2 = 0
        
        # Converter para numérico
        allele1_num = allele1 + 0
        allele2_num = allele2 + 0
        
        total++
        depth_sum += $15
        
        # Calcular diferença entre alelos
        diff = (allele2_num > allele1_num) ? allele2_num - allele1_num : allele1_num - allele2_num
        diff_sum += diff
        
        # Contar homozigotos vs heterozigotos
        if (int(allele1_num + 0.5) == int(allele2_num + 0.5)) {
            hom++
        } else {
            het++
        }
        
        # Imprimir linha com NaN substituído
        print sample, $1, $2, $3, $4, allele1_num, allele2_num, $7, $8, $9, $10, $11, $12, $13, $14, $15
    }
    END {
        if (total > 0) {
            printf "%s\t%d\t%.2f\t%.2f\t%d\t%d\n", sample, total, depth_sum/total, diff_sum/total, hom, het
        }
    }
    ' "$TSV_FILE")
    
    # Extract statistics line
    stats_line=$(echo "$awk_result" | tail -1)
    variant_data=$(echo "$awk_result" | head -n -1)
    
    # Add to consolidated file (uncompressed)
    if [ -n "$variant_data" ]; then
        echo "$variant_data" >> "$output_file"
    fi
    
    # Add statistics to summary
    if [ -n "$stats_line" ] && echo "$stats_line" | grep -q $'\t'; then
        echo "$stats_line" >> "$sample_summary"
        sample_variants=$(echo "$stats_line" | cut -f2)
        total_variants=$((total_variants + sample_variants))
    fi
    
    variant_count_num=$(echo "$variant_data" | wc -l 2>/dev/null)
    if [ -z "$variant_count_num" ]; then
        variant_count_num=0
    fi
    echo "OK - $variant_count_num variants"
done

echo "======================================"
echo "CONSOLIDATION COMPLETE"
echo "======================================"
echo "OK - Output file: $output_file"
echo "OK - Sample summary: $sample_summary"
echo ""
echo "SUMMARY STATISTICS:"
echo "  Total samples processed: $counter"
echo "  Total variants: $total_variants"
if [ $counter -gt 0 ]; then
    avg_variants=$((total_variants / counter))
    echo "  Average variants per sample: $avg_variants"
fi

echo ""
echo "FILE SIZE:"
if [ -f "$output_file" ]; then
    output_size=$(du -h "$output_file" | cut -f1)
    echo "  Consolidated TSV: $output_size ($(wc -l < "$output_file") lines)"
fi
if [ -f "$sample_summary" ]; then
    summary_size=$(du -h "$sample_summary" | cut -f1)
    echo "  Sample summary: $summary_size ($(wc -l < "$sample_summary") lines)"
fi
