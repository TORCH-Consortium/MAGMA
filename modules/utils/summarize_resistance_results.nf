process UTILS_SUMMARIZE_RESISTANCE_RESULTS {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("minor_variants/*")
        path("major_variants/*")
        path("structural_variants/*")

    output:
        path("resistance_summaries")

    script:
       
        """
        summarize_resistance.py major_variants minor_variants structural_variants resistance_summaries
        """

    stub: 

        """
        mkdir resistance_summaries
        """ 

}
