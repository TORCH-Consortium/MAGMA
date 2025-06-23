process FASTQ_VALIDATOR {
    tag "${sampleName}"
    // tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(sampleName), path(sampleRead)
    val ready


    // tuple val(meta), path(sampleRead)



    output:
    path("*.fastq_report.csv")                                       , emit: fastq_report
    tuple val(sampleName), path(sampleRead)                          , emit: reads

    // tuple val(meta), path("*.bam")   , emit: bam
    // path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fastqvalidator \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqvalidator: \$(fastqvalidator --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${sampleRead.simpleName}.check.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqvalidator: \$(fastqvalidator --version)
    END_VERSIONS
    """
}
