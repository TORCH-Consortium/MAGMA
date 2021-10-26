process DELLY_CALL {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bai), path(recalibratedBam)
    path(reference)

    output:
    tuple val(sampleName), path("*.delly.bcf")

    script:

    """
    ${params.delly_path} call \\
        -g ${reference} \\
        ${recalibratedBam} \\
        ${params.arguments} \\
        -o ${sampleName}.delly.bcf

    """

    stub:

    """
    echo "${params.delly_path} call \\
        -g ${reference} \\
        ${recalibratedBam} \\
        ${params.arguments} \\
        -o ${sampleName}.delly.bcf"

    touch ${sampleName}.delly.bcf
    """

}
