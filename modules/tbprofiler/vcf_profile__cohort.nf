process TBPROFILER_VCF_PROFILE__COHORT {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(mergedVcfIndex), path(mergedVcf)
    path(resistanceDb)

    // output:


    script:
    def optionalDb  = resistanceDb ? "--db ${resistanceDb}" : ""

    """
    ${params.tbprofiler_path} vcf_profile  \\
        ${optionalDb} \\
        ${params.arguments} \\
        ${mergedVcf}

    """

    stub:
    """
    """

}
