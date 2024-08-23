process UTILS_REFORMAT_LOFREQ {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    container "ghcr.io/torch-consortium/magma/misc:2.0.0-alpha"

    input:
        tuple val(sampleName), path(lofreqVcf)

    output:
        tuple val(sampleName), path("*lofreq.reformat.corrected.vcf")

    script:
       
        """
        reformat_lofreq.py ${lofreqVcf} \\
            ${sampleName} \\
            ${sampleName}.lofreq.reformat.vcf

        reduce_strand_bias.py \\
            ${params.cutoff_strand_bias} \\
            ${sampleName}.lofreq.reformat.vcf  \\
            ${sampleName}.lofreq.reformat.corrected.vcf 
        """

    stub: 

        """
        touch ${sampleName}.lofreq.reformat.vcf
        touch ${sampleName}.lofreq.reformat.corrected.vcf 
        """ 

}
