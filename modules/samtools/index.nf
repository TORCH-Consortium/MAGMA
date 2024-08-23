process SAMTOOLS_INDEX {
    tag "${sampleName}"
    label 'cpu_2_memory_2'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bam)

    output:
        tuple val(sampleName), path("*.bai"), path(bam)

    script:

        """
        ${params.samtools_path} index ${bam}
        """

    stub:

        """
        echo "samtools index ${bam}"
        touch ${sampleName}.bam
        touch ${sampleName}.bam.bai
        """

}
