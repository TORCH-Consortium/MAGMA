process LOFREQ_CALL {
    tag "${sampleName}"
    label 'cpu_8_memory_4'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bai), path(dindleBam)
        path(ref_fasta)
        path("*")

    output:
        tuple val(sampleName), path("*.LoFreq.vcf")

    script:

        """
        ${params.lofreq_path} call-parallel \\
            -f ${ref_fasta} \\
            --pp-threads ${task.cpus} \\
            ${params.arguments} \\
            ${dindleBam} \\
            -o ${sampleName}.LoFreq.vcf

        # Trigger the process again by chaging this script
        """

    stub:

        """
        echo "lofreq call \\
            -f ${ref_fasta} \\
            ${params.arguments} \\
            --call-indels \\
            ${dindleBam} \\
            > ${sampleName}.LoFreq.vcf"

        touch ${sampleName}.LoFreq.vcf
        """

}
