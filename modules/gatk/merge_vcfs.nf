nextflow.enable.dsl = 2


params.results_dir = "${params.outdir}/gatk4/merge_vcfs"
params.save_mode = 'copy'
params.should_publish = true

process GATK_MERGE_VCFS {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    gatk MergeVcfs -Xmx${task.memory.giga}G \\
    -I $JOINT_NAME/$JOINT_NAME.filtered_SNP_exc-rRNA.vcf.gz \\
    -I $JOINT_NAME/$JOINT_NAME.raw_INDEL.vcf.gz \\
    -O $JOINT_NAME/$JOINT_NAME.filtered_SNP.RawIndels.vcf.gz
    """

    stub:

    """

    """
}

