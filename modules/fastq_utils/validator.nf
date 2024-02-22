process FASTQ_VALIDATOR {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)

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
        else
            VALIDATED=0
            STATUS="passed"
        fi

        echo -e "!{sampleName}\t${VALIDATED}" >> !{sampleName}.check.${STATUS}.tsv

        '''

    stub: 

        """
        touch ${sampleName}.check.tsv 
        """ 

}
