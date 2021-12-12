process UTILS_QUANTTB_SAMPLE_QC {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish, pattern: "*quanttb.qc.csv"

    input:
        tuple val(sampleName), path(quanttbStats)
        val(relabundanceThreshold)
        tuple val(sampleName), val(bamRgString), path(sampleReads)

    output:
        path("*quanttb.qc.csv") ,emit: quanttb_sample_qc_ch
        tuple val(sampleName), val(bamRgString), path("*quanttb.qc.csv"), path(sampleReads) ,emit: qc_samplereads_tuple

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
