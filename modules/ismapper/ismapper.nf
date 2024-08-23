process ISMAPPER {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        tuple val(sampleName), val(meta), path(reads)
        path(ref_gbk)
        path(ref_fasta)
        path(query_multifasta)


    output:
        tuple val(sampleName), path("*.ismapper.vcf"), emit: formatted_vcf
        tuple val(sampleName), path("ismapper"), emit: results_dir

    script:

        """

        ismap \\
            --t $task.cpus \\
            --output_dir $sampleName \\
            --queries $query_multifasta \\
            --reference $ref_gbk \\
            --reads $reads

        # Reorganize output files
        mkdir ismapper
        mv $sampleName/*/* ismapper/


        # Reformat the is-mapper output to a standard VCF
        # FIXME the IS element name is hard-coded
        convert_ismapper_to_vcf.py \\
            --is_mapper_dir ismapper/IS6110/ \\
            --reference_file $ref_fasta \\
            --query_file $query_multifasta \\
            --output_vcf_file ${sampleName}.ismapper.vcf

        """
}
