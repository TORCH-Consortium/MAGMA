process NTMPROFILER_PROFILE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish, pattern: "vcf/*vcf.gz", saveAs: { f -> f.tokenize("/").last()}


    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)

    //output:
        //tuple val(sampleName), path("vcf/*"), emit: per_sample_results

    script:

        """
        ${params.ntmprofiler_path} profile \\
            -1 ${sampleReads[0]} \\
            -2 ${sampleReads[1]} \\
            -p ${sampleName}
        """

    stub:
        """
        touch ${sampleName}.txt
        """
}


