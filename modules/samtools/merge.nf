process SAMTOOLS_MERGE {
    tag "${sampleName}"
    label 'cpu_medium_memory_medium'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path("bams/*")

    output:
        tuple val(sampleName), path("*.sorted_reads.bam")

    script:

        """
        ${params.samtools_path} merge \\
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

        touch ${sampleName}.sorted_reads.bam
        """

}
