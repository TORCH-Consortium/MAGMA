nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/collect_wgs_metrics"
params.save_mode = 'copy'
params.should_publish = true
params.arguments = " --READ_LENGTH 0 --COVERAGE_CAP 10000 --COUNT_UNPAIRED"

process GATK_COLLECT_WGS_METRICS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bam)
    path(ref_fasta)

    output:
    tuple val(sampleName), path(".*WgsMetrics.txt")


    script:

    """
    gatk CollectWgsMetrics -Xmx${task.memory.giga}G \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        ${arguments} \\
        -O ${sampleName}.WgsMetrics.txt
    """

    stub:

    """
    echo "gatk CollectWgsMetrics -Xmx${task.memory.giga}G \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        ${arguments} \\
        -O ${sampleName}.WgsMetrics.txt"

    touch ${sampleName}.WgsMetrics.txt
    """
}

