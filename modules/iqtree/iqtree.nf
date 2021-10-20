process IQTREE {
    tag "${joint_name}"
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(incComplexFasta)

    output:
    //FIXME

    script:

    """
    ${params.iqtree_path} -s ${incComplexFasta} \\
        -T ${task.cpus} \\
        ${params.arguments} \\
        --prefix ${joint_name}.95X.IncComplex
    """

    stub:

    """
    """

}
