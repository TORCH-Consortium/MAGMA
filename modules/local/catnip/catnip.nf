process CATNIP {

    tag "${prefix} - ${snp_threshold} SNP"

    publishDir params.results_dir,
        mode: params.save_mode,
        enabled: params.should_publish

    input:
    path(snp_matrix)
    val(snp_threshold)
    val(prefix)

    output:
    path("${prefix}.${snp_threshold}SNPcluster.tsv"),
        emit: cluster_annotation

    script:
    """
    python3 catnip.py \
        ${snp_matrix} \
        ${snp_threshold} \
        > ${prefix}.${snp_threshold}SNPcluster.tsv
    """

    stub:
    """
    touch ${prefix}.${snp_threshold}SNPcluster.tsv
    """
}
