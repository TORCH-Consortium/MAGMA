/*
FIXME: Documentation comments

*/



process BGZIP {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(annotatedVcf)

    output:
    tuple val(joint_name), path("*.gz")

    script:

    """
    bgzip ${annotatedVcf}
    """

    stub:

    """
    touch ${joint_name}.annotated.vcf.gz
    """

}