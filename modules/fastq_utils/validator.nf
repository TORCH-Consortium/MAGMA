process FASTQ_VALIDATOR {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)

    output:
        tuple val(sampleName), path("*.check.tsv") 

    shell:
       
        '''
        !{params.fastq_validator_path} !{sampleReads} \\
            2>!{sampleName}.command.log || true

        cp !{sampleName}.command.log .command.log


        TEMP=$(tail -n 1 !{sampleName}.command.log)

        if [ "$(echo "$TEMP")" == "OK" ]; then
            VALIDATED=1
        else
            VALIDATED=0
        fi

        echo -e "!{sampleName}\t!{bamRgString}\t${VALIDATED}" > !{sampleName}.check.tsv

        '''

    stub: 

        """
        touch ${sampleName}.check.tsv 
        """ 

}
