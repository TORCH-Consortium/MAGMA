process UTILS_SUMMARIZE_RESISTANCE_RESULTS {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("xbs_output_dir/smartt/analyses/drug_resistance/minor_variants/*")
        path("xbs_output_dir/smartt/analyses/drug_resistance/major_variants/*")
        path(lit_csv) 
        path(lit_forced_csv) 

    output:
        path("resistance_summaries")

    script:
       
        """
        summarize_resistance.py xbs_output_dir smartt lit.csv lit_forced.csv resistance_summaries
        """

    stub: 

        """
        mkdir resistance_summaries
        """ 

}
