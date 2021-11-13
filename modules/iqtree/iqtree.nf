process IQTREE {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    val(prefix)
    tuple val(joint_name), path(fasta)

    output:
    tuple val(joint_name),
        path("joint.${prefix}.bionj"),
        path("joint.${prefix}.ckp.gz"),
        path("joint.${prefix}.iqtree"),
        path("joint.${prefix}.log"),
        path("joint.${prefix}.mldist"),
        path("joint.${prefix}.model.gz"),
        path("joint.${prefix}.treefile")

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
    touch joint.${prefix}.bionj
    touch joint.${prefix}.ckp.gz
    touch joint.${prefix}.iqtree
    touch joint.${prefix}.log
    touch joint.${prefix}.mldist
    touch joint.${prefix}.model.gz
    touch joint.${prefix}.treefil
    """

}
