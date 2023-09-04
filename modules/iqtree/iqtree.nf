process IQTREE {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

        val(prefix)
        tuple val(joint_name), path(fasta)

    output:

        tuple val(joint_name), path("${joint_name}.${prefix}.treefile"), emit: tree_tuple
        path("${joint_name}.${prefix}.bionj")
        path("${joint_name}.${prefix}.ckp.gz")
        path("${joint_name}.${prefix}.iqtree")
        path("${joint_name}.${prefix}.log")
        path("${joint_name}.${prefix}.mldist")
        path("${joint_name}.${prefix}.model.gz")
        path("${joint_name}.${prefix}.treefile")


    script:
        if(params.iqtree_standard_bootstrap) {
            arguments = '-bb 1000'
        } else if(params.iqtree_fast_ml_only) {
            arguments = '-fast'
        } else if(params.iqtree_fast_bootstrapped_phylogeny) {
            arguments = '-bb 1000 -alrt 1000'
        } else if(params.iqtree_accurate_ml_only) {
            arguments = '-allnni'
        } else {
        //NOTE: Use iqtree_accurate_ml_only as the default
            arguments = '-allnni'
        }

        """
        ${params.iqtree_path} \\
            -s ${fasta} \\
            -T AUTO \\
            ${arguments} \\
            --prefix ${joint_name}.${prefix}
        """

    stub:

        """
        touch ${joint_name}.${prefix}.bionj
        touch ${joint_name}.${prefix}.ckp.gz
        touch ${joint_name}.${prefix}.iqtree
        touch ${joint_name}.${prefix}.log
        touch ${joint_name}.${prefix}.mldist
        touch ${joint_name}.${prefix}.model.gz
        touch ${joint_name}.${prefix}.treefil
        """

}
