process GATK_COMBINE_GVCFS {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        val(joint_name)
        val(gvcfs_string)
        path(gvcfs)
        path(ref_fasta)
        path(ref_exit_rif_gvcf)
        path(ref_exit_rif_gvcf_tbi) //NOTE: The optional refExitRifGvcfTbi had to be provided as a separate input since it was causing comparison error when staged via  path("*")
        path("*")

    output:
        tuple val(joint_name),  path("*.combined.vcf.gz.tbi"), path("*.combined.vcf.gz")


    script:

        def optionalRefExitRifGvcf  = ( params.use_ref_gvcf ) ? " --variant ${ref_exit_rif_gvcf} " : ""

        """
        ${params.gatk_path} CombineGVCFs --java-options "-Xmx${task.memory.giga}G" \\
            -R ${ref_fasta} \\
            ${params.arguments} \\
            --variant ${gvcfs_string} \\
            ${optionalRefExitRifGvcf} \\
            -O ${joint_name}.combined.vcf.gz

            cp ${joint_name}.combined.vcf.gz .command.out
        """

    stub:

        """
        touch ${joint_name}.combined.vcf.gz
        touch ${joint_name}.combined.vcf.gz.tbi
        """
}

