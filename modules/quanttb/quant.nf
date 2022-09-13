process QUANTTB_QUANT {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish, pattern: "output/*txt"

    input:
        tuple val(sampleName), val(bamRgString) ,path(sampleReads)

    output:
        tuple val(sampleName), path("output/*.quant.txt")         ,emit: quanttb_report_tuple
        tuple val(sampleName), val(bamRgString), path(sampleReads)  ,emit: samplereads_tuple

    script:

        """
        sleep 10
        ${params.quanttb_path} quant -f ${sampleReads} -o ${sampleName}.quant.txt -k
        sleep 10
        """

    stub:

        """

        echo "quanttb quant -f ${sampleReads} -o ${sampleName}.quant.txt -k"

        mkdir output
        touch output/${sampleName}.quant.txt
        """

}
