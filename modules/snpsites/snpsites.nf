process SNPSITES {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    val(prefix)
    tuple val(joint_name), path(alignmentFasta)

    output:
    tuple val(joint_name), path("*.variable.${prefix}.fa")


    script:

    """
    ${params.snpsites_path} -o ${joint_name}.variable.${prefix}.fa ${alignmentFasta}
    """

    stub:

    """
    touch ${joint_name}.variable.${prefix}.fa
    """

}
