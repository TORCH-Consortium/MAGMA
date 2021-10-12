nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/samtools/index"
params.save_mode = 'copy'
params.should_publish = true


process SAMTOOLS_MERGE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bam)

    output:
    tuple val(sampleName), path("*.bai")

    script:

    """
    samtools index ${bam}
    """

    stub:

    """
    echo "samtools index ${bam}"
    """

}
