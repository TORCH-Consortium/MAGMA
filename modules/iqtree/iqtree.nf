process IQTREE {
    tag "${joint_name}"
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(fasta)

    // output:
    //FIXME

    script:

    """
    ${params.iqtree_path} \\
        -s ${fasta} \\
        -T ${task.cpus} \\
        ${params.arguments} \\
        --prefix ${joint_name}.${prefix}
    """

    stub:

    """
    """

}
