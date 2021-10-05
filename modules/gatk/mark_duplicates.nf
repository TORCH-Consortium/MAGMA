nextflow.enable.dsl = 2

params.gatk_path = "gatk"
params.java_opts = "-Xms4000m"
params.compression_level = 5

process GATK_MARK_DUPLICATES {
    tag "${sampleId}"
    label 'gatk4_container'

    input:
    val(sampleId)
    path(input_mapped_merged_bam)

    output:
    val(sampleId)
    path("${sampleId}_merged.deduped.bam")
    path("${sampleId}_merged.deduped.metrics.txt")

    script:

    """
    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opts}" \
                        MarkDuplicates \
                        --INPUT ${input_mapped_merged_bam} \
                        --OUTPUT ${sampleId}_merged.deduped.bam \
                        --METRICS_FILE ${sampleId}_merged.deduped.metrics.txt \
                        --VALIDATION_STRINGENCY SILENT \
                        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                        --ASSUME_SORT_ORDER "queryname" \
                        --CREATE_MD5_FILE true
    """

    stub:

    """
    touch ${sampleId}_merged.deduped.bam
    touch ${sampleId}_merged.deduped.metrics.txt
    """

}
