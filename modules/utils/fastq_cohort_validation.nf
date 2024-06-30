process UTILS_FASTQ_COHORT_VALIDATION {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("fastq_reports/*")
        path(magma_validated_samplesheet_json)

    output:
        path("magma_analysis.json"), emit: magma_analysis_json

    script:

        """
        csvtk concat fastq_reports/* |  csvtk csv2json -k file > merged_fastq_reports.json

        fastq_cohort_validation.py ${magma_validated_samplesheet_json} merged_fastq_reports.json magma_analysis.json

        rm merged_fastq_reports.json

        """


    stub:

        """
        touch ${params.vcf_name}.passed.tsv ${params.vcf_name}.failed.tsv
        """

}
