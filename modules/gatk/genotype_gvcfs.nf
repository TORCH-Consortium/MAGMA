/*
FIXME: Documentation comments

*/


process GATK_GENOTYPE_GVCFS {
    tag "${joint_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    tuple val(joint_name), path(combinedVcf)
    path(ref_fasta)


    output:
    tuple val(joint_name), path("*.raw_variants.vcf.gz")


    script:

    """
    gatk GenotypeGVCFs -Xmx${task.memory.giga}G \\
        -R ${ref_fasta} \\
        -V ${combinedVcf} \\
        ${params.arguments} \\
        -O ${joint_name}.raw_variants.vcf.gz
    """

    stub:

    """
    """
}

