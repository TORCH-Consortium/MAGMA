process SNP_SITES {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(alignmentFasta)
    val(type)


    output:
    tuple val(joint_name), path("*.95X.variable.${type}.fa")


    script:

    """
    ${params.snp_sites_path} -o ${joint_name}.95X.variable.${type}.fa ${alignmentFasta}
    """

    stub:

    """
    touch ${joint_name}.95X.variable.${type}.fa
    """

}
