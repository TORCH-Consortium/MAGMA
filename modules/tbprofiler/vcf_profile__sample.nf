process TBPROFILER_VCF_PROFILE__SAMPLE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(resistanceDb)
    tuple val(sampleName), path(lofreqVcf)

    output:
    //FIXME Find the ideal output


    script:
    """
    tb-profiler vcf_profile \\
        --lofreq_sample_name ${sampleName} \\
        --db ${resistanceDb} \\
        ${lofreqVcf}
    """

    stub:
    """
    """

}
