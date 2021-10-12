nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/base_recalibrator"
params.save_mode = 'copy'
params.should_publish = true

process GATK_BASE_RECALIBRATOR {
    tag "$sampleName"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(dedupedBam)
    path(dbsnp)
    path(ref_fasta)

    output:
    tuple val(sampleName), path(".*recal_data.table"), path(dedupedBam)

    script:

    """
    gatk BaseRecalibrator -Xmx${task.memory.giga}G \\
        --known-sites ${dbsnp} \\
        -R ${ref_fasta} \\
	    -I ${dedupedBam} \\
        -O ${sampleName}.recal_data.table
    """

    stub:

    """
    echo "gatk BaseRecalibrator -Xmx${task.memory.giga}G \\
        --known-sites ${dbsnp} \\
        -R ${ref_fasta} \\
	    -I ${dedupedBam} \\
        -O ${sampleName}.recal_data.table"

    touch ${sampleName}.recal_data.table
    """
}
