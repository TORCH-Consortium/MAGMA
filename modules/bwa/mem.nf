process BWA_MEM {
    tag "${sampleName}"
    label 'cpu_8_memory_8'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), val(meta), path(sampleReads)
        path(reference)
        path("*")

    output:
        tuple val(sampleName), path("*.sorted_reads.bam")


    script:

        """
        ${params.bwa_path} mem \\
            -M \\
            ${params.arguments} \\
            -t ${task.cpus} \\
            -R "${meta.bam_rg_string}"  \\
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
            -R ${meta.bam_rg_string}  \\
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
