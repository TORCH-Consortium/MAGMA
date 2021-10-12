nextflow.enable.dsl = 2


params.results_dir = "${params.outdir}/fastqc"
params.save_mode = 'copy'
params.should_publish = true

process FASTQC {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(sampleReads)

    output:
    tuple path('*.html'), path('*.zip')


    script:

    """
    fastqc \\
        ${sampleReads} \\
        -t ${task.cpus}
    """

    stub:
    """
    echo "fastqc \\
              ${sampleReads} \\
              -t ${task.cpus}"

    touch ${sampleName}.html

    touch ${sampleName}.zip
    """
}
