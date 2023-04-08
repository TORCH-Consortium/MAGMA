process TBPROFILER_PROFILE__BAM {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path("*.bai"), path(bam)
        path(resistanceDb)

    output:
        tuple val(sampleName), path("vcf/*"), emit: per_sample_results

    script:
        def optionalDb  = resistanceDb ? "--db ${resistanceDb}" : ""

        """
        ${params.tbprofiler_path} profile \\
            -a ${bam} \\
            ${optionalDb} \\
            -p ${sampleName}
        """

    stub:
        """
        touch ${sampleName}.txt
        """
}


