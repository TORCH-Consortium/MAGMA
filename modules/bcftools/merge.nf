process BCFTOOLS_MERGE {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(vcfs_string_ch)
        path("*")

    output:
        tuple val(params.vcf_name), path("*.${params.file_format}.vcf")

    script:

        """
        bcftools merge -o ${params.vcf_name}.${params.file_format}.vcf ${vcfs_string_ch}
        """

    stub:

        """
        touch ${params.vcf_name}.LoFreq.vcf
        """

}
