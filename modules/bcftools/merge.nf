process BCFTOOLS_MERGE {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(vcfs_string_ch)
        path("*")

    output:
        tuple val(params.vcf_name), path("*.${params.file_format}.vcf.gz"), path("*.vcf.gz.csi")


    script:

        """
        bcftools merge -o ${params.vcf_name}.${params.file_format}.vcf ${vcfs_string_ch}
        bgzip ${params.vcf_name}.${params.file_format}.vcf
        ${params.bcftools_path} index ${params.vcf_name}.${params.file_format}.vcf.gz
        """

    stub:

        """
        touch ${params.vcf_name}.${params.file_format}.vcf.gz
        touch ${params.vcf_name}.${params.file_format}.vcf.gz.tbi

        """

}
