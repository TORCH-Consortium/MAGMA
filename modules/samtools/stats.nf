process SAMTOOLS_STATS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bam)
    path(reference)

    output:
    tuple val(sampleName), path(".*SamtoolStats.txt")

    script:

    """
    ${params.samtools_path} stats \\
        ${params.arguments} \\
        ${bam} \\
        -r ${reference} \\
    > ${sampleName}.SamtoolStats.txt
    """

    stub:

    """
    touch ${sampleName}.SamtoolStats.txt
    """

}
