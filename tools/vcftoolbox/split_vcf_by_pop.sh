#!/bin/bash

#Exit on error
set -e

vcf_input="$1"
vcf_names="$2"
indpop_file="$3"

vcf_dir="vcf_split_directory"

##### Verify that the output directory exists #####
if [[ ! -d "${vcf_dir}" ]]; then
    echo "ERROR: Failed to create output VCF directory: ${vcf_dir}" >&2
    exit 1
fi

##### Check input files #####
if [[ -z "$vcf_input" ]]; then
    echo "ERROR: VCF file is not provided." >&2
    exit 1
fi

# Verify that input VCF contains at least one variant
if ! bcftools view -H "$vcf_input" | head -n 1 | grep -q .; then
    echo "ERROR: Input VCF contains no variant records."
    exit 1
fi

if [[ -z "$indpop_file" ]]; then
    echo "ERROR: Population file is empty or not provided." >&2
    exit 1
fi

##### Check indpop_file structure #####
local line_num=0
local error_found=0
while IFS= read -r line || [[ -n "$line" ]]; do
    ((line_num++))

    # Check empty lines
    [[ -z "$line" ]] && continue

    # Ensure exactly two tab-separated columns (one tab delimiter only)
    local col_count
    col_count=$(echo "$line" | awk -F'\t' '{print NF}')
    if [[ "$col_count" -ne 2 ]]; then
        echo "ERROR: Population files does not contain exactly 2 tab-separated columns" >&2
        error_found=1
    fi

done < "$indpop_file"

if [[ "$error_found" -eq 1 ]]; then
    echo "ERROR: Population file has an invalid structure. Expected: 2 tab-separated columns and no header." >&2
    exit 1
fi

echo "INFO: indpop_file structure OK."


##### Global variables ######
pop_file_list="${vcf_dir}/population_files_list.txt" #Temporary file in working directory

#########################################################
# Function: split_individuals_by_pop
# Description: Creates lists of individuals by population
#########################################################

split_individuals_by_pop(){
    local indpop_file="$1"
    
    ##### Check if file exists #####
    if [[ ! -f "$indpop_file" ]]; then
        echo "Error: population file not found: $indpop_file" >&2
        exit 1
    fi

    ##### Store individuals grouped by population #####
    #Identify unique populations
    declare -A pop_inds

    while IFS=$'\t' read -r ind pop; do
        #Skip empty lines
        [[ -z "$ind" || -z "$pop" ]] && continue

        #Initialize array if not exists
        if [[ -z "${pop_inds[$pop]}" ]]; then
            pop_inds[$pop]="" ## Initialize string that will store individuals for this population
        fi

        #Store individual
        pop_inds[$pop]+="${ind}"$'\n'

    done < "$indpop_file"

    ##### Initialize the list file (empty the file contents if it already exists) #####
    > "$pop_file_list"

    ##### Write one individual list file per population #####
    for pop in "${!pop_inds[@]}"; do
                local output_file="${vcf_dir}/Ind_list_${pop}.txt"
                echo -n "${pop_inds[$pop]}" > "$output_file"
                echo "$output_file|$pop" >> "$pop_file_list"
    done

    ##### Number of populations detected #####
    local pop_count="${#pop_inds[@]}"
    echo "INFO: ${pop_count} population(s) detected: ${!pop_inds[*]}"
    
    ###### Clear associative array ######
    unset pop_inds
    
}

#########################################################
# Function: split_vcf_by_pop
# Description: Creates VCF file for each population
#########################################################

split_vcf_by_pop() {
    ##### Parameters #####
    local vcf="$1"
    local original_name="$2"

    #### Check if VCF file exists #####
    if [[ ! -f "$vcf" ]]; then
        echo "Error: VCF file not found: $vcf" >&2
        exit 1
    fi

    ##### Check if population files list exists #####
    if [[ ! -f "$pop_file_list" ]]; then
        echo "Error: Population files list not found: $pop_file_list" >&2
        exit 1
    fi

    ##### Counter for created VCF files #####
    local vcf_created=0

    ##### Process each individual list file #####
    while IFS='|' read -r ind_list pop_name; do
        [[ ! -f "$ind_list" ]] && continue
        
        # Count individuals in the list
        local ind_count=$(grep -c ^ "$ind_list")
        
        # Extract base name
        local vcf_basename
        regex='\(([^)]+)\)[[:space:]]*$'
        if [[ "$original_name" =~ $regex ]]; then
            #Extract content between last parentheses
            vcf_basename="${BASH_REMATCH[1]}"
        else
            # No parentheses, use original name
            vcf_basename=$(basename "$original_name")
        fi
    
        #Remove file extension
        vcf_basename=${vcf_basename%.*}
        
        local output_vcf="${vcf_dir}/${vcf_basename}_${pop_name}.vcf"
        
        # Create a temporary gzipped VCF for indexing
        local tmp_vcf="${output_vcf}.gz"
        
        # Subset VCF and compress
        bcftools view -S "$ind_list" --force-samples "$vcf" -Oz -o "$tmp_vcf"
        
        # Verify that VCF creation succeeded
        if [[ $? -eq 0 && -f "$tmp_vcf" ]]; then
            
            # Index the output VCF
            bcftools index -t "$tmp_vcf"
            
            # Convert compressed VCF back to uncompressed VCF format
            bcftools view "$tmp_vcf" -Ov -o "$output_vcf"
            
            # Cleanup temporary files
            rm -f "$tmp_vcf" "${tmp_vcf}.tbi"
            
            ((vcf_created ++))
        else
            echo "Failed to create VCF for $pop_name"
        fi
        
        
        done < "$pop_file_list"

        ##### Verify that filtered VCF is not empty ######
        if [[ ! -s "$output_vcf" ]]; then
            echo "Output VCF not created: $output_vcf" >&2
            exit 1
        fi

        if ! bcftools view -H "$output_vcf" | head -n 1 | grep -q .; then
            echo "ERROR: Filtered VCF contains no variants."
            exit 1
        fi   

        echo "INFO: ${vcf_created} VCF file(s) successfully created."
        
}

########################################
# Main execution
########################################
main(){
    local vcf_input="$1"
    local vcf_names="$2"
    local indpop_file="$3"

    split_individuals_by_pop "$indpop_file"
    split_vcf_by_pop "$vcf_input" "$vcf_names"
}

##### Execution #####
main "$vcf_input" "$vcf_names" "$indpop_file"

#Cleanup temporary files
rm -f "${vcf_dir}"/Ind_list_*.txt
rm -f "${pop_file_list}"