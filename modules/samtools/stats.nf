nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/samtools/stats"
params.save_mode = 'copy'
params.should_publish = true

params.arguments = "-F DUP,SUPPLEMENTARY,SECONDARY,UNMAP,QCFAIL"

process SAMTOOLS_STATS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bam)
    path(ref_fasta)

    output:
    tuple val(sampleName), path(".*SamtoolStats.txt")

    script:

    """
    samtools stats \\
        ${arguments} \\
        ${bam} \\
        -r ${ref_fasta} \\
    > ${sampleName}.SamtoolStats.txt
    """

    stub:

    """
    touch ${sampleName}.SamtoolStats.txt
    """

}
