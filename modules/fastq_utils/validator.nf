process FASTQ_VALIDATOR {
    tag "${magmaName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    stageInMode 'copy'
    maxRetries 3
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }


    input:
        tuple val(magmaName), path(sampleRead)
        val ready

    output:
        tuple val(magmaName), path("*.check.*tsv")
        path("*.check.*tsv")                          , emit: check_result
        tuple val(magmaName), path(sampleRead)        , emit: reads

    shell:

        '''
        seqkit stats -a -T  !{sampleRead}  > !{sampleRead}.seqkit_out.txt
        cat *seqkit_out.txt | csvtk space2tab | csvtk tab2csv > !{sampleRead}.seqkit_stats.csv

        md5sum !{sampleRead} > !{sampleRead}.md5sum_out.txt
        cat *md5sum_out.txt | csvtk space2tab | csvtk tab2csv | csvtk add-header -n md5sum,file > !{sampleRead}.md5sum_stats.csv

        du -shL !{sampleRead} > !{sampleRead}.du_out.txt
        cat *du_out.txt | csvtk tab2csv | csvtk add-header -n size,file > !{sampleRead}.du_stats.csv



        csvtk join -f file \\
        !{sampleRead}.seqkit_stats.csv \\
        !{sampleRead}.md5sum_stats.csv \\
        !{sampleRead}.du_stats.csv \\
        > !{sampleRead}.fastq_statistics.csv


        rm *_out.txt *_stats.csv


        !{params.fastq_validator_path} !{sampleRead} \\
        2>!{sampleRead}.command.log || true

        cp !{sampleRead}.command.log .command.log



        TEMP=$(tail -n 1 !{sampleRead}.command.log)


        if [ "$(echo "$TEMP")" == "OK" ]; then
            VALIDATED=1
            STATUS="passed"
            echo -e "!{sampleRead}\t${VALIDATED}" > !{sampleRead}.check.${STATUS}.tsv
            exit 0

        else
            VALIDATED=0
            STATUS="failed"
            echo -e "!{sampleRead}\t${VALIDATED}" > !{sampleRead}.check.${STATUS}.tsv
            exit 1
        fi


        '''

    stub:

        """
        touch ${sampleRead}.check.tsv
        """

}
