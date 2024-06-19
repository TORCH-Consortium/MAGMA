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
        seqkit stats -a -T  *fastq*  > ${sampleName}.seqkit_stats.csv
        cat ${sampleName}.seqkit_stats.csv | csvtk space2tab | csvtk tab2csv > ${sampleName}.seqkit_stats.final.csv

        md5sum *fastq* > ${sampleName}.md5sum_stats.csv 
        cat ${sampleName}.md5sum_stats.csv | csvtk space2tab | csvtk tab2csv | csvtk add-header -n md5sum,file > ${sampleName}.md5sum_stats.final.csv

        du -shL *fastq* > ${sampleName}.du_stats.csv 
        cat ${sampleName}.du_stats.csv | csvtk tab2csv | csvtk add-header -n size,file > ${sampleName}.du_stats.final.csv

        """

    stub: 

        """
        touch ${sampleName}.check.tsv 
        """ 

}
