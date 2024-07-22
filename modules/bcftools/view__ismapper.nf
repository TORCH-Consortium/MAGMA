process BCFTOOLS_VIEW__ISMAPPER {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(vcf)

    output:
        tuple val(sampleName), path("*.ismapper.vcf.gz"), path("*.ismapper.vcf.gz.csi")


    script:

        """
        bcftools view -Oz -o ${sampleName}.ismapper.vcf.gz $vcf

        bcftools index ${sampleName}.ismapper.vcf.gz
        """

    stub:

        """
        touch ${sampleName}.ismapper.vcf.gz
        """

}
