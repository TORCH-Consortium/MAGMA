process UTILS_QUANTTB_SAMPLE_QC {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish, pattern: "*quanttb.qc.csv"

    input:
        tuple val(sampleName), path(quanttbStats)
        val(relabundanceThreshold)

    output:
        tuple val(sampleName), path("*quanttb.qc.csv")


    script:
        """
        quanttb_qc.py \\
            --input ${quanttbStats} \\
            --output ${sampleName}.quanttb.qc.csv \\
            --relative_abundance_threshold ${relabundanceThreshold} \\
            --write_header false \\
            --derived_sample_name ${sampleName}
        """
}
