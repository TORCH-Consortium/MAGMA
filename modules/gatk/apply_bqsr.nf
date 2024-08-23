process GATK_APPLY_BQSR {
    tag "$sampleName"
    label 'cpu_2_memory_2'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(recalibrationTable), path(dedupedBam)
        path(ref_fasta)
        path("*")

    output:
        tuple val(sampleName), path("*.baserecalibrated.bam")

    script:

        """
        ${params.gatk_path} ApplyBQSR --java-options "-Xmx${task.memory.giga}G" \\
            -R ${ref_fasta} \\
            -I ${dedupedBam} \\
            --bqsr ${recalibrationTable} \\
            -O ${sampleName}.baserecalibrated.bam
        """

    stub:

        """
        echo "gatk ApplyBQSR -Xmx${task.memory.giga}G \\
            -R ${ref_fasta} \\
            -I ${dedupedBam} \\
            --bqsr ${recalibrationTable} \\
            -O ${sampleName}.baserecalibrated.bam"

        touch ${sampleName}.baserecalibrated.bam
        """
}
