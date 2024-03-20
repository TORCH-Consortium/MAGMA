process NTMPROFILER_PROFILE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish, pattern: "vcf/*vcf.gz", saveAs: { f -> f.tokenize("/").last()}


    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)
        path(resistanceDb)

    //output:
        //tuple val(sampleName), path("vcf/*"), emit: per_sample_results

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


