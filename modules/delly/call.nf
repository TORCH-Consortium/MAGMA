process DELLY_CALL {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(recalibratedBam)
    path(ref_fasta)

    output:
    tuple val(sampleName), path("*.delly.bcf")

    script:

    """
    delly call \\
        -g ${ref_fasta} \\
        ${recalibratedBam} \\
        ${arguments} \\
        -o ${sampleName}.delly.bcf

    """

    stub:

    """
    echo "delly call \\
        -g ${ref_fasta} \\
        ${recalibratedBam} \\
        ${arguments} \\
        -o ${sampleName}.delly.bcf"

    touch ${sampleName}.delly.bcf
    """

}
