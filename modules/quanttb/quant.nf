process QUANTTB_QUANT {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    beforeScript "source ${params.quanttb_venv}/activate.sh"
    //NOTE: Uncomment if there is a corresponding activate script which you'd like to run
    // afterScript = "source ${params.quanttb_venv}/deactivate.sh"
    afterScript "deactivate"

    input:
    tuple val(sampleName), path(sampleReads)

    output:
    tuple val(sampleName), path("*.quant.txt")

    script:

    """
    ${params.quanttb_path} quant ${sampleReads} -o ${sampleName}.quant.txt -k
    """

    stub:

    """

    echo "quanttb quant ${sampleReads} -o ${sampleName}.quant.txt -k"

    touch ${sampleName}.quant.txt
    """

}
