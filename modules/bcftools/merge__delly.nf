process BCFTOOLS_MERGE__DELLY {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path vcf_files

    output:
        tuple val(params.vcf_name), path("*.vcf.gz.csi"), path("*.${params.file_format}.vcf.gz")

    script:
        """
        # Extract sample prefixes
        prefixes=(\$(for file in ${vcf_files.join(' ')}; do basename \$file | cut -d '.' -f 1; done | sort -u))

        # Concatenate files for each sample prefix
        concat_files=()
        for prefix in \${prefixes[@]}; do
            bcftools concat \${prefix}.*.vcf.gz -o \${prefix}.concat.vcf
            bgzip \${prefix}.concat.vcf
            ${params.bcftools_path} index \${prefix}.concat.vcf.gz
            concat_files+=("\${prefix}.concat.vcf.gz")
        done

        # Merge the concatenated files
        bcftools merge \${concat_files[@]} -o ${params.vcf_name}.${params.file_format}.vcf
        bgzip ${params.vcf_name}.${params.file_format}.vcf
        ${params.bcftools_path} index ${params.vcf_name}.${params.file_format}.vcf.gz
        """

    stub:
        """
        touch ${params.vcf_name}.${params.file_format}.vcf.gz
        touch ${params.vcf_name}.${params.file_format}.vcf.gz.csi
        """
}
