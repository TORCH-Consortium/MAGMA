process UTILS_REFORMAT_LOFREQ {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(lofreqVcf)

    output:
        tuple val(sampleName), path("*LoFreq.Reformat.vcf")

    script:
       
        """
        reformat_lofreq.py ${lofreqVcf} \\
            ${sampleName} \\
            ${sampleName}.LoFreq.Reformat.vcf
        """

    stub: 

        """
        touch ${sampleName}.LoFreq.Reformat.vcf
        """ 

}
