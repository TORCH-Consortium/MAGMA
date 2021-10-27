process GATK_COLLECT_WGS_METRICS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bam)
    path(reference)

    output:
    tuple val(sampleName), path("*.WgsMetrics.txt")


    script:

    """
    ${params.gatk_path} CollectWgsMetrics --java-options "-Xmx${task.memory.giga}G" \\
        -R ${reference} \\
        -I ${bam} \\
        ${params.arguments} \\
        -O ${sampleName}.WgsMetrics.txt
    """

    stub:

    """
    echo "gatk CollectWgsMetrics -Xmx${task.memory.giga}G \\
        -R ${reference} \\
        -I ${bam} \\
        ${params.arguments} \\
        -O ${sampleName}.WgsMetrics.txt"

    touch ${sampleName}.WgsMetrics.txt
    """
}

