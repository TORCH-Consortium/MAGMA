process BWA_MEM {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), val(bamRgString), path(sampleReads)
    path(reference)
    path("*")

    output:
    tuple val(sampleName), path("*.sorted_reads.bam")


    script:

    """
    ${params.bwa_path} mem \\
        -M \\
        -t ${task.cpus} \\
        -R "${RG}"  \\
        ${reference} \\
        ${sampleReads} \\
    | ${params.samtools_path} sort \\
        -@ ${task.cpus} \\
        -O BAM \\
        -o ${sampleName}.sorted_reads.bam -
    """

    stub:

    """
    echo "${params.bwa_path} mem \\
        -M \\
        -t ${task.cpus} \\
        -R ${RG}  \\
        ${reference} \\
        ${sampleReads} \\
    | ${params.samtools_path} sort \\
        -@ ${task.cpus} \\
        -O BAM \\
        -o ${sampleName}.sorted_reads.bam -
    "

    touch ${sampleName}.sorted_reads.bam
    """

}
