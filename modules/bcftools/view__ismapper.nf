process BCFTOOLS_VIEW__ISMAPPER {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*.ismapper.vcf.gz")

    script:

        """
        bcftools view -Oz -o ${sampleName}.ismapper.vcf.gz $vcf

        """

    stub:

        """
        touch ${sampleName}.ismapper.vcf.gz
        """

}
