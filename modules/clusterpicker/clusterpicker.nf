nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/cluster_picker"
params.save_mode = 'copy'
params.should_publish = true


process CLUSTER_PICKER{
    tag ""
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(fasta)
    path(newick_tree)
    val(bootstrap_1)
    val(bootstrap_2)
    val(max_cluster_size)
    val(algorithm)

    output:

    script:
    //NOTE: The java_opts are propogated via bash as per the wrapper script in bioconda
    // https://github.com/bioconda/bioconda-recipes/blob/master/recipes/clusterpicker/cluster-picker.sh


    //NOTE: As per the manual of cluster-picker
    // cluster-picker input-fasta.fas input-tree.nwk bootstrap bootstrap genetic-distance maxClusterSize dist


    // cluster-picker $OUT_DIR/fasta/$JOINT_NAME/$JOINT_NAME.95X.variable.IncComplex.5SNPcluster.fa $OUT_DIR/phylogeny/$JOINT_NAME/$JOINT_NAME.95X.IncComplex.5SNPcluster.treefile 0 0 $FRACTION1  0 gap

    """
    jvm_mem_opts="-Xmx${task.memory.giga}G"

    cluster-picker \\
        ${fasta} \\
        ${newick_tree} \\
        ${bootstrap_1} \\
        ${bootstrap_2} \\
        ${max_cluster_size} \\
        ${algorithm}
    """

    stub:

    """
    """

}
