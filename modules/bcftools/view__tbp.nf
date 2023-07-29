process BCFTOOLS_VIEW__TBP {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bcf)

    output:
        tuple val(sampleName), path("*.bcf.gz"), path("*.bcf.gz.csi")

    shell:

        '''
        !{params.bcftools_path} view !{params.arguments}  !{sampleName}.delly.bcf -o !{sampleName}.filtered.delly.bcf
        bgzip !{sampleName}.filtered.delly.bcf
        !{params.bcftools_path} index !{sampleName}.filtered.delly.bcf.gz
        '''

    stub:

        """
        touch ${sampleName}.bcf.gz
        touch ${sampleName}.bcf.gz.csi
        """

}
