process CATNIP {

    tag "${joint_name} - ${snp_threshold} SNP"

    publishDir params.results_dir,
        mode: params.save_mode,
        enabled: params.should_publish

    input:
    tuple val(joint_name), path(snp_matrix)
    val(snp_threshold)
    val(prefix)

    output:
    tuple val(joint_name),
          path("${joint_name}.${prefix}.${snp_threshold}SNPcluster.tsv"),
          emit: cluster_annotation

    script:
    """
    python ${projectDir}/bin/catnip.py \
        ${snp_matrix} \
        ${snp_threshold} \
        > ${joint_name}.${prefix}.${snp_threshold}SNPcluster.tsv
    """

    stub:
    """
    touch ${joint_name}.${prefix}.${snp_threshold}SNPcluster.tsv
    """
}
