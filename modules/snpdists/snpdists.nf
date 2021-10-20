process SNPDISTS {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(alignmentFasta)
    val(type)

    output:
    path("*.snp_dists.tsv")

    script:

    """
    ${params.snp_dists_path} ${alignmentFasta} -b \\
    > ${joint_name}.${type}.snp_dists.tsv
    """

    stub:

    """
    touch ${joint_name}.${type}.snp_dists.tsv
    """

}
