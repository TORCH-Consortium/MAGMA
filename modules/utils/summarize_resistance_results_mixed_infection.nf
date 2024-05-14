process UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("minor_variants/*")
        path("structural_variants/*")

    output:
        path("combined_resistance_summaries_mixed_infection_samples")

    script:
       
        """
        summarize_resistance_mixed_infection.py minor_variants structural_variants combined_resistance_summaries_mixed_infection_samples
        """

    stub: 

        """
        mkdir combined_resistance_summaries_mixed_infection
        """ 

}
