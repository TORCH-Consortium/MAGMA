process CATNIP {

    tag "${prefix} - ${snp_threshold} SNP"

    publishDir params.results_dir,
        mode: params.save_mode,
        enabled: params.should_publish

    input:
    tuple val(joint_name), path(snp_matrix), path(treefile)
    val(snp_threshold)
    val(prefix)
    path(catnip_script)

    output:
    tuple val(joint_name),
          path("${joint_name}.${prefix}.${snp_threshold}SNPcluster.tsv"),
          emit: cluster_annotation
    
    tuple val(joint_name),
          path("${joint_name}.${prefix}.${snp_threshold}SNPcluster.nexus"),
          emit: nexus_tree

    tuple val(joint_name),
          path("*.cluster.txt"),
          emit: sample_cluster_files

    script:
    """
    python3 ${catnip_script} \
        ${snp_matrix} \
        ${joint_name}.${prefix}.${snp_threshold}SNPcluster.tsv \
        --threshold ${snp_threshold} \
        --tree ${treefile} \
        --tree-out ${joint_name}.${prefix}.${snp_threshold}SNPcluster.nexus
    """

    stub:
    """
    touch ${prefix}.${snp_threshold}SNPcluster.tsv
    touch ${joint_name}.${prefix}.${snp_threshold}SNPcluster.nexus
    """
}
