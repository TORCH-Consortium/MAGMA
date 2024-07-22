process BCFTOOLS_MERGE__LOFREQ {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path(vcfs_file)
        path("*")

    output:
        tuple val(params.vcf_name), path("*.vcf.gz.csi"), path("*.${params.file_format}.vcf.gz")

    script:

        """
        bcftools merge -o ${params.vcf_name}.${params.file_format}.vcf -l ${vcfs_file}
        bgzip ${params.vcf_name}.${params.file_format}.vcf
        ${params.bcftools_path} index ${params.vcf_name}.${params.file_format}.vcf.gz
        """

    stub:

        """
        touch ${params.vcf_name}.${params.file_format}.vcf.gz
        touch ${params.vcf_name}.${params.file_format}.vcf.gz.csi
        """

}
