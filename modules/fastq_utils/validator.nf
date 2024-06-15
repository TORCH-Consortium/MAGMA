process FASTQ_VALIDATOR {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    stageInMode 'copy'
    maxRetries 3
    //errorStrategy 'retry'
    //NOTE: Default action is to ignore the process if the second attempt fails
    errorStrategy { (task.attempt <= task.maxRetries) ? 'retry' : 'ignore' }



    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)
        val ready

    output:
        tuple val(sampleName), path("*.check.*tsv") 
        path("*.check.*tsv"), emit: check_result

    shell:
       
        '''
        !{params.fastq_validator_path} !{sampleReads} \\
            2>!{sampleName}.command.log || true

        cp !{sampleName}.command.log .command.log


        TEMP=$(tail -n 1 !{sampleName}.command.log)

        if [ "$(echo "$TEMP")" == "OK" ]; then
            VALIDATED=1
            STATUS="passed"
            echo -e "!{sampleName}\t${VALIDATED}" > !{sampleName}.check.${STATUS}.tsv
            exit 0
        else
            VALIDATED=0
            STATUS="failed"
            echo -e "!{sampleName}\t${VALIDATED}" > !{sampleName}.check.${STATUS}.tsv
            exit 1
        fi


        '''

    stub: 

        """
        touch ${sampleName}.check.tsv 
        """ 

}
