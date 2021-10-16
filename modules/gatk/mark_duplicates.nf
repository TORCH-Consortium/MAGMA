nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/mark_duplicates"
params.save_mode = 'copy'
params.should_publish = true


process GATK_MARK_DUPLICATES {
    tag "$sampleName"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(mergedBam)

    output:
    tuple val(sampleName), path(".*dedup_reads.bam")
    tuple val(sampleName), path(".*MarkDupMetrics.txt")

    script:

    """
    gatk MarkDuplicates -Xmx${task.memory.giga}G \\
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
