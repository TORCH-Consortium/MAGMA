process TBPROFILER_VCF_PROFILE__LOFREQ {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(lofreqVcf)
        path(resistanceDb)

    output:
        tuple val(sampleName), path("results/*")
        path("results/*"), emit: resistance_json


    script:
        def optionalDb  = resistanceDb ? "--db ${resistanceDb.name}" : ""

//FIXME
        """
        ${params.tbprofiler_path} vcf_profile \\
            ${optionalDb} \\
            ${mergedLofreqVcf}
        """

    stub:
        """
        mkdir results
        touch results/${sampleName}.results.json
        """

}
