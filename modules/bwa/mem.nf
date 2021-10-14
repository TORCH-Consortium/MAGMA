nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/bwa/mem"
params.save_mode = 'copy'
params.should_publish = true

process BWA_MEM {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(sampleReads), val(RG)
    path(reference)

    output:
    path("*.sorted_reads.bam")


    script:

    """
    bwa mem \\
        -M \\
        -t ${task.cpus} \\
        -R ${RG}  \\
        ${reference} \\
        ${sampleReads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -O BAM \\
        -o ${sampleName}.sorted_reads.bam -
    """

    stub:

    """
    echo "bwa mem \\
        -M \\
        -t ${task.cpus} \\
        -R ${RG}  \\
        ${reference} \\
        ${sampleReads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -O BAM \\
        -o ${sampleName}.sorted_reads.bam -
    "

    touch ${sampleName}.sorted_reads.bam
    """

}
