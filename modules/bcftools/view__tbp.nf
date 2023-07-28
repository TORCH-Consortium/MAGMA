process BCFTOOLS_VIEW__TBP {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bcf)

    //output:
    //    tuple val(sampleName), path("*.potentialSV.vcf.gz")

    shell:

        '''
        !{params.bcftools_path} view !{params.arguments}  !{sampleName}.delly.bcf -o !{sampleName}_filtered.delly.bcf
        bgzip !{sampleName}_filtered.delly.bcf
        !{params.bcftools_path} index !{sampleName}_filtered.delly.bcf.gz
        '''

    stub:

        """
        touch ${sampleName}.potentialSV.vcf.gz
        """

}
