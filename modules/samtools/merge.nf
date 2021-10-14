nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/samtools/merge"
params.save_mode = 'copy'
params.should_publish = true


process SAMTOOLS_MERGE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path("bams/*")

    output:
    tuple val(sampleName), path(".*sorted_reads.bam")

    script:

    """
    samtools merge \\
        -f \\
        ${sampleName}.sorted_reads.bam \\
        bams/* \\
        -@ ${task.cpus}
    """

    stub:

    """
    echo "samtools merge \\
        -f \\
        ${sampleName}.sorted_reads.bam \\
        bams/* \\
        -@ ${task.cpus}"

    """

}
