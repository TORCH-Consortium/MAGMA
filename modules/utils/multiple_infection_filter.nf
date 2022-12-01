process UTILS_MULTIPLE_INFECTION_FILTER {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("collated_results/*")

    output:
        path("approved_samples.relabundance.tsv"), emit: approved_samples
        path("rejected_samples.relabundance.tsv"), emit: rejected_samples

    script:
       
        """
        multiple_infection_filter.py --indir collated_results --relative_abundance_threshold ${params.rel_abundance_cutoff}
        """

    stub: 

        """
        touch approved_samples.relabundance.tsv
        touch rejected_samples.relabundance.tsv
        """ 

}
