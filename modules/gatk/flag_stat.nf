process GATK_FLAG_STAT {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bam)
    path(ref_fasta)
    path("*")


    output:
    tuple val(sampleName), path("*.FlagStat.txt")

    script:

    """
    ${params.gatk_path} FlagStat --java-options "-Xmx${task.memory.giga}G" \\
        -R ${ref_fasta} \\
        -I ${bam}  \\
    > ${sampleName}.FlagStat.txt

    """

    stub:

    """
    echo "gatk FlagStat -Xmx${task.memory.giga}G \\
        -R ${ref_fasta} \\
        -I ${bam}  \\
    > ${sampleName}.FlagStat.txt"


    touch ${sampleName}.FlagStat.txt
    """
}

