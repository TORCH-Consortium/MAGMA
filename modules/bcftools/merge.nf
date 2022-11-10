process BCFTOOLS_MERGE {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bcf)

    output:
        tuple val(sampleName), path("*.FIXME")

    script:

        """

        """

    stub:

        """
        touch ${sampleName}.potentialSV.vcf.gz
        """

}
