process LOFREQ_CALL {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bai), path(dindleBam)
    path(ref_fasta)
    path("*")

    output:
    tuple val(sampleName), path("*.LoFreq.vcf")

    script:

    """
    ${params.lofreq_path} call \\
        -f ${ref_fasta} \\
        ${params.arguments} \\
        ${dindleBam} \\
    > ${sampleName}.LoFreq.vcf
    """

    stub:

    """
    echo "lofreq call \\
        -f ${ref_fasta} \\
        ${params.arguments} \\
        --call-indels \\
        ${dindleBam} \\
        > ${sampleName}.lofreq.vcf"

    touch ${sampleName}.lofreq.vcf
    """

}
