process UTILS_FASTQ_STATS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    container "ghcr.io/torch-consortium/magma/misc:2.0.0-alpha"

    input:
        tuple val(sampleName), path(sampleReads)

//    output:
//        tuple val(sampleName), path("*.check.*tsv"), emit: check_result

    script:
       
        """
        seqkit stats -a -T  *fastq*  > ${sampleName}.seqkit.stats.csv

        md5sum *fastq* > ${sampleName}.md5sum.stats.csv

        du -sh *fastq* > ${sampleName}.du.stats.csv

        """

    stub: 

        """
        touch ${sampleName}.check.tsv 
        """ 

}
