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
        trap "rm -rf \$tmp_dir" EXIT

        # Extract sample prefixes
        prefixes=(\$(for file in ${vcf_files.join(' ')}; do basename \$file | cut -d '.' -f 2; done | sort -u))

        # Concatenate files for each sample prefix
        concat_files=()
        for prefix in \${prefixes[@]}; do
            files=(\$(ls \${prefix}.*.bcf.gz \${prefix}.*.vcf.gz 2>/dev/null))
            if [ \${#files[@]} -gt 0 ]; then
                concat_file="\${tmp_dir}/\${prefix}.concat.vcf"
                bcftools concat \${files[@]} -o \${concat_file}
                bgzip \${concat_file}
                sorted_file="\${tmp_dir}/\${prefix}.sorted.vcf.gz"
                bcftools sort \${concat_file}.gz -o \${sorted_file}
                ${params.bcftools_path} index \${sorted_file}
                concat_files+=("\${sorted_file}")
            fi
        done

        # Merge the concatenated files
        if [ \${#concat_files[@]} -gt 0 ]; then
            bcftools merge \${concat_files[@]} -o ${params.vcf_name}.${params.file_format}.vcf
            bgzip ${params.vcf_name}.${params.file_format}.vcf
            ${params.bcftools_path} index ${params.vcf_name}.${params.file_format}.vcf.gz
        fi
        """

    stub:
        """
        touch ${params.vcf_name}.${params.file_format}.vcf.gz
        touch ${params.vcf_name}.${params.file_format}.vcf.gz.csi
        """
}
