process UTILS_SUMMARIZE_RESISTANCE_RESULTS {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("minor_variants/*")
        path("major_variants/*")
        path("resources/lit.csv") 
        path("resources/lit_forced.csv") 

    output:
        path("resistance_summaries")

    script:
       
        """
        summarize_resistance.py path("") smartt_master_branch path("resistance_summaries")
        """

    stub: 

        """
        mkdir resistance_summaries
        """ 

}
