process GATK_MARK_DUPLICATES {
    tag "$sampleName"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(mergedBam)

    output:
    tuple val(sampleName), path(".*dedup_reads.bam"),              emit: bam_tuple
    tuple val(sampleName), path(".*MarkDupMetrics.txt"),           emit: metrics

    script:

    """
    ${params.gatk_path} MarkDuplicates -Xmx${task.memory.giga}G \\
        --METRICS_FILE ${sampleName}.MarkDupMetrics.txt \\
        -I ${mergedBam} \\
        -O ${sampleName}.dedup_reads.bam
    """

    stub:

    """
    echo " gatk MarkDuplicates -Xmx${task.memory.giga}G \\
        --METRICS_FILE ${sampleName}.MarkDupMetrics.txt \\
        -I ${mergedBam} \\
        -O ${sampleName}.dedup_reads.bam"

    touch ${sampleName}.MarkDupMetrics.txt
    touch ${sampleName}.dedup_reads.bam
    """

}
