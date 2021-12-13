process CLUSTERPICKER {
    tag "snpCount - ${snpCount}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(joint_name), path(fasta), path(newickTree)
        val(snpCount)
        val(prefix)

    output:
        tuple val(joint_name), path("*${snpCount}SNPcluster*")

    shell:
        //NOTE: The java_opts are propogated via bash as per the wrapper script in bioconda
        // https://github.com/bioconda/bioconda-recipes/blob/master/recipes/clusterpicker/cluster-picker.sh


        //NOTE: As per the instructions in the cluster-picker manual
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

        cp joint.!{prefix}_clusterPicks_log.txt joint.!{prefix}_!{snpCount}SNPclusterPicks_log.txt
        cp joint.!{prefix}_clusterPicks.nwk joint.!{prefix}_!{snpCount}SNPclusterPicks.nwk
        cp joint.!{prefix}_clusterPicks.nwk.figTree joint.!{prefix}_!{snpCount}SNPclusterPicks.nwk.figTree
        cp joint.variable.!{prefix}.fa_joint.!{prefix}_clusterPicks.fas joint.variable.!{prefix}.fa_joint.!{prefix}_!{snpCount}SNPclusterPicks.fas

        '''

    stub:

        """
        """

}
