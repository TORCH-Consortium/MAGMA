process UTILS_FASTQ_STATS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)
        val ready

    output:
        tuple val(sampleName), path("*.check.*tsv") 
        path("*.check.*tsv"), emit: check_result

    shell:
       
        '''
        seqkit stats -a -T  !{sampleReads} | csvtk pretty -t  > ${sampleName}.seqkit.stats.csv

        md5sum !{} >> 

        du -h !{}

        '''

    stub: 

        """
        touch ${sampleName}.check.tsv 
        """ 

}
