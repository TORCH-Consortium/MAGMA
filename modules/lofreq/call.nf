/*
FIXME: Documentation comments

*/


process LOFREQ_CALL {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(dindleBam), path(bai)
    path(ref_fasta)

    output:
    tuple val(sampleName), path("*.lofreq.vcf")

    script:

    """
    lofreq call \\
        ${arguments} \\
        --call-indels \\
        ${dindleBam} \\
    > ${sampleName}.lofreq.vcf
    """

    stub:

    """
    echo "lofreq call \\
        ${arguments} \\
        --call-indels \\
        ${dindleBam} \\
        > ${sampleName}.lofreq.vcf"

    touch ${sampleName}.lofreq.vcf
    """

}
