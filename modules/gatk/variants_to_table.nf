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
        !{params.gatk_path} VariantsToTable --java-options "-Xmx!{task.memory.giga}G" \\
            -V !{vcf} \\
            !{params.arguments} \\
            -O /dev/stdout \\
        | sed -e 's/^\\t//g' \\
        | sed -e 's/*/-/g' \\
        | sed -e 's/\\./-/g' \\
        | sed '2,${/^.*\\(-.*\\)\\{'"!{params.median_coverage_cutoff}"',\\}.*$/d}' \\
        | !{params.datamash_path} transpose \\
        | sed -e 's/^/>/g' \\
        | sed -e 's/-GT/\\n/g' \\
        | sed -e 's/\\t//g' \\
        > !{joint_name}.!{prefix}.fa


        '''

    stub:

        """
        touch ${joint_name}.${prefix}.fa
        """
}

