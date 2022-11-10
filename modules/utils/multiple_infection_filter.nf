process UTILS_MULTIPLE_INFECTION_FILTER {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        indir
        relative_abundance_threshold

    output:
        path("approved_samples.relabundance.tsv")
        path("rejected_samples.relabundance.tsv")

    script:
       
        """
        multiple_infection_filter.py         
        """

    stub: 

        """
        touch 
        """ 

}
