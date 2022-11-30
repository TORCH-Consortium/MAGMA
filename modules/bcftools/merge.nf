process BCFTOOLS_MERGE {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("*")

    output:
        path("*.LoFreq.vcf.gz")

    script:

        """
        bcftools merge -o ${sampleName}.LoFreq.vcf.gz
        """

    stub:

        """
        touch ${sampleName}.LoFreq.vcf.gz
        """

}
