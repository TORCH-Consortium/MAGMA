process ISMAPPER {

    input:
        tuple val(sampleName), val(meta), path(reads)
        path(reference)
        path(query)


    output:
        tuple val(meta), path("ismapper/*"), emit: results

    script:
        query_base = query.getSimpleName()

        """

        ismap \\
            --t $task.cpus \\
            --output_dir $sampleName \\
            --queries $query \\
            --reference $reference \\
            --reads $reads

        # Reorganize output files
        mkdir ismapper
        mv $sampleName/*/* ismapper/


        # FIXME This would need to updated later for accommodating various static lengths per insertion element. Work with a CSV file name,element_sequence,sequence_length
        #convert_ismapper_to_vcf.py

        """
}
