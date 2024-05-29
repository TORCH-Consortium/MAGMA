process UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path(merge_cohort_stats)
        path("minor_variants/*")
        path("structural_variants/*")

    output:
        path("combined_resistance_summaries_mixed_infection_samples"), optional: true

    script:
       
        """
        summarize_resistance_mixed_infection.py ${merge_cohort_stats} minor_variants structural_variants combined_resistance_summaries_mixed_infection_samples
        """

    stub: 

        """
        mkdir combined_resistance_summaries_mixed_infection
        """ 

}
