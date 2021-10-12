nextflow.enable.dsl = 2


//NOTE: FASTQC
// - The SEQ_R1 and SEQ_R2 are processed separately to accomodate each others absence (XBS_map#L26)
// - Integrity check is done using `gzip -t`


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
    fastqc  -t ${task.cpus}
    """

    stub:
    """
    touch ${sampleName}.html

    touch ${sampleName}.zip
    """
}
