process ISMAPPER {

    input:
        tuple val(sampleName), val(meta), path(reads)
        path(ref_gbk)
        path(ref_fasta)
        path(query_multifasta)


    output:
        tuple val(meta), path("ismapper/*"), emit: results

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


        # FIXME This would need to updated later for accommodating various static lengths per insertion element. Work with a CSV file name,element_sequence,sequence_length
        convert_ismapper_to_vcf.py \\
            --is_mapper_file \\
            --reference_file $ref_fasta \\
            --query_file $query_multifasta \\
            --output_vcf_file output.vcf

        """
}
