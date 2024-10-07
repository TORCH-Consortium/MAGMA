process UTILS_MERGE_STRUCTURAL_VARIANTS {
    //publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("fastq_reports/*")
        path(magma_validated_samplesheet_json)

    script:

        """
        merge_structural_vcfs.py \\
            --delly_vcf \\
            --ismapper_vcf \\
            --bed_file \\
            --existing_json_file \\
            --cleaned_json_file \\
            --output_vcf

        """


}
