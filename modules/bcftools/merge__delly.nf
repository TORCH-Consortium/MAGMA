process BCFTOOLS_MERGE__DELLY {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    stageInMode 'copy'

    input:
        path vcf_files

    output:
        tuple val(params.vcf_name), path("*.vcf.gz.csi"), path("*.${params.file_format}.vcf.gz")

    script:
        """
        bcftools_merge__delly.sh
        """

    stub:
        """
        touch ${params.vcf_name}.${params.file_format}.vcf.gz
        touch ${params.vcf_name}.${params.file_format}.vcf.gz.csi
        """
}
