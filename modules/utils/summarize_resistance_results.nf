process UTILS_SUMMARIZE_RESISTANCE_RESULTS {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("minor_variants/*")
        path("major_variants/*")
        path(lit_csv) 
        path(lit_forced_csv) 

    output:
        path("resistance_summaries")

    script:
       
        """
        summarize_resistance.py  
        """

    stub: 

        """
        mkdir resistance_summaries
        """ 

}
