process CLUSTERPICKER {
    tag "snpCount - ${snpCount}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(fasta), path(newickTree)
    val(snpCount)

    output:
    //FIXME Find out what's the ideal output

    shell:
    //NOTE: The java_opts are propogated via bash as per the wrapper script in bioconda
    // https://github.com/bioconda/bioconda-recipes/blob/master/recipes/clusterpicker/cluster-picker.sh


    //NOTE: As per the manual of cluster-picker
    // cluster-picker input-fasta.fas                                    input-tree.nwk                                  bootstrap bootstrap genetic-distance max_cluster_size dist
    // cluster-picker $JOINT_NAME.95X.variable.IncComplex.5SNPcluster.fa $JOINT_NAME.95X.IncComplex.5SNPcluster.treefile 0         0         $FRACTION1       0              gap

    '''

    cp !{fasta} !{fasta.getBaseName()}.!{snpCount}SNPcluster.fa

    cp !{newickTree} !{newickTree.getBaseName()}.!{snpCount}SNPcluster.treefile

    # computes the genetic-distance argument
    FRACTION=$(echo !{snpCount}/$(cat !{fasta} | head -n 2 | tail -n 1| wc -c) | bc -l)

    jvm_mem_opts="-Xmx!{task.memory.giga}G"

    !{params.clusterpicker_path} \\
        !{fasta} \\
        !{newickTree} \\
        !{params.bootstrap_1} \\
        !{params.bootstrap_2} \\
        $FRACTION \\
        !{params.max_cluster_size} \\
        !{params.algorithm}
    '''

    stub:

    """
    """

}
