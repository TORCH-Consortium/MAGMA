#!/usr/bin/env bash

# Create temporary directory
tmp_dir="interim"

mkdir $tmp_dir

# Extract unique sample prefixes from filenames
prefixes=( $(for file in *.bcf.gz *.vcf.gz; do basename "$file" | cut -d '.' -f 2; done | sort -u) )

# Process each sample prefix
for prefix in "${prefixes[@]}"; do
    # Find related files for the prefix
    files=( "*${prefix}.[^.]*.bcf.gz" "*${prefix}.[^.]*.bcf.gz" )

    echo $files

    if [[ ${#files[@]} -gt 0 ]]; then
        # Concatenate files
        concat_file="${tmp_dir}/${prefix}.concat.vcf.gz"
        bcftools concat ${files[@]} -o "$concat_file"
        bcftools view "$concat_file" -Oz -o "$concat_file"

        # Sort and index the concatenated file
        sorted_file="${tmp_dir}/${prefix}.concat.sorted.vcf.gz"
        bcftools sort "$concat_file" -Oz -o "$sorted_file"
        bcftools index "$sorted_file"

        # Add sorted file to merged list
        concat_files+=("$sorted_file")
    fi
done

# Merge concatenated files (if any)
if [[ ${#concat_files[@]} -gt 0 ]]; then
    bcftools merge "${concat_files[@]}" -o joint.delly.vcf
    bgzip joint.delly.vcf
    bcftools index joint.delly.vcf.gz
fi

# NOTE For now commenting this for debugging.
# Clean up temporary directory
#rm -rf "$tmp_dir"
