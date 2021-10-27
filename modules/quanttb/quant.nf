process QUANTTB_QUANT {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(sampleReads)

    output:
    tuple val(sampleName), path("output/*.quant.txt")

    script:

    """
    ${params.quanttb_path} quant -f ${sampleReads} -o ${sampleName}.quant.txt -k
    """

    stub:

    """

    echo "quanttb quant -f ${sampleReads} -o ${sampleName}.quant.txt -k"

    touch ${sampleName}.quant.txt
    """

}
