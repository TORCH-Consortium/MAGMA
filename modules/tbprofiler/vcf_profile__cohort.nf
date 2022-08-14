process TBPROFILER_VCF_PROFILE__COHORT {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(joint_name), path(mergedVcfIndex), path(mergedVcf)
        path(resistanceDb)

    output:
        path("results/*")


    script:
        def optionalDb  = ( ref_exit_rif_gvcf.simpleName != "NONE") ? "--db ${resistanceDb.name}" : ""

        def optionallyLoadLibraryForContainers = (optionalDb != "") ? "cd ${resistanceDb}; ${params.tbprofiler_path} load_library ${resistanceDb.name}; cd ../" : ""

        """
        ${optionallyLoadLibraryForContainers}

        ${params.tbprofiler_path} vcf_profile  \\
            ${optionalDb} \\
            ${mergedVcf}

        """

    stub:
        """
        """

}
