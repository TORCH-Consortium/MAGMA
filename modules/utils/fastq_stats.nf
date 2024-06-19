process UTILS_FASTQ_STATS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    container "ghcr.io/torch-consortium/magma/misc:2.0.0-alpha"

    input:
        tuple val(sampleName), path(sampleReads)

    output:
        path("*fastq_stats.csv")

    script:
       
        """
        seqkit stats -a -T  *fastq*  > ${sampleName}.seqkit.txt
        cat ${sampleName}.seqkit.txt | csvtk space2tab | csvtk tab2csv > ${sampleName}.seqkit_stats.final.csv

        md5sum *fastq* > ${sampleName}.md5sum.txt 
        cat ${sampleName}.md5sum.txt | csvtk space2tab | csvtk tab2csv | csvtk add-header -n md5sum,file > ${sampleName}.md5sum_stats.csv

        du -shL *fastq* > ${sampleName}.du.txt 
        cat ${sampleName}.du.txt | csvtk tab2csv | csvtk add-header -n size,file > ${sampleName}.du_stats.csv



        csvtk join -f file \\
            ${sampleName}.seqkit_stats.final.csv \\
            ${sampleName}.md5sum_stats.csv \\
            ${sampleName}.du_stats.csv \\
        > ${sampleName}.fastq_stats.csv


        """

    stub: 

        """
        touch ${sampleName}.check.tsv 
        """ 

}
