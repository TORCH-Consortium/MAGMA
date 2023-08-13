process GATK_VARIANTS_TO_TABLE {
    tag "${joint_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        val(prefix)
        tuple val(joint_name), path(vcfIndex), path(vcf)

    output:
        tuple val(joint_name), path("*.fa")


    shell:

        '''
        !{params.gatk_path} VariantsToTable --java-options "-Xmx!{task.memory.giga}G" \
            -V !{vcf} !{params.arguments} \
            -O !{joint_name}.!{prefix}.table \
        | variant_table_to_fasta.py !{joint_name}.!{prefix}.table \
            !{joint_name}.!{prefix}.fa \
            !{params.site_representation_cutoff}
        '''

    stub:

        """
        touch ${joint_name}.${prefix}.fa
        """
}

