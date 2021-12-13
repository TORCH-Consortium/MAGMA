process SNPDISTS {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(prefix)
        tuple val(joint_name), path(alignmentFasta)

    output:
        tuple val(joint_name), path("*.snp_dists.tsv")

    script:

        """
        ${params.snpdists_path} ${alignmentFasta} -b \\
        > ${joint_name}.${prefix}.snp_dists.tsv
        """

    stub:

        """
        touch ${joint_name}.${prefix}.snp_dists.tsv
        """

}
