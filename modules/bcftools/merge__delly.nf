/*
 * Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
 *
 * This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
 *
 * For quick overview of GPL-3 license, please refer
 * https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
 *
 * - You MUST keep this license with original authors in your copy
 * - You MUST acknowledge the original source of this software
 * - You MUST state significant changes made to the original software
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program . If not, see <http://www.gnu.org/licenses/>.
 */
process BCFTOOLS_MERGE__DELLY {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path vcf_files

    output:
        tuple val(params.vcf_name), path("*.vcf.gz.csi"), path("*.${params.file_format}.vcf.gz")

    script:
        """

        # Create temporary directory
        tmp_dir=\$(mktemp -d -t bcftools-XXXXXX)

        # Extract unique sample prefixes from filenames
        prefixes=( \$(for file in *.bcf.gz *.vcf.gz; do basename $file | cut -d '.' -f 2; done | sort -u) )

        # Process each sample prefix
        for prefix in "${prefixes[@]}"; do
        # Find related files for the prefix
        files=( "${prefix}"*.bcf.gz "${prefix}"*.vcf.gz )

        if [[ ${#files[@]} -gt 0 ]]; then
            # Concatenate files
            concat_file="${tmp_dir}/${prefix}.concat.vcf"
            bcftools concat "${files[@]}" -o "$concat_file"
            bgzip "$concat_file"

            # Sort and index the concatenated file
            sorted_file="${tmp_dir}/${prefix}.sorted.vcf.gz"
            bcftools sort "$concat_file.gz" -o "$sorted_file"
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

        # Clean up temporary directory
        rm -rf "$tmp_dir"

        """

    stub:
        """
        touch ${params.vcf_name}.${params.file_format}.vcf.gz
        touch ${params.vcf_name}.${params.file_format}.vcf.gz.csi
        """
}
