process NTMPROFILER_PROFILE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish, pattern: "vcf/*vcf.gz", saveAs: { f -> f.tokenize("/").last()}


    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)

    output:
        path("results"),  emit: per_sample_results
        path("results/*json"), emit: profile_json

    script:

        """
        ${params.ntmprofiler_path} profile \\
            -1 ${sampleReads[0]} \\
            -2 ${sampleReads[1]} \\
            -p ${sampleName}    \\
            -d results    \\
            --txt 
        """

    stub:
        """
        mkdir ${sampleName}
        touch "${sampleName}/${sampleName}.json"
        """
}


