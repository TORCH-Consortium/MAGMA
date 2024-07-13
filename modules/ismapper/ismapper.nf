process ISMAPPER {

    input:
        tuple val(sampleName), val(meta), path(reads)
        path(reference)
        path(query)


    output:
        tuple val(meta), path("ismapper/*"), emit: results

    script:
        def reference_name = reference.getName().replace(".gz", "")
        def query_name = query.getName().replace(".gz", "")
        query_base = query.getSimpleName()

        """
        if [ "$ref_compressed" == "true" ]; then
            gzip -c -d $reference > $reference_name
        fi
        if [ "$query_compressed" == "true" ]; then
            gzip -c -d $query > $query_name
        fi

        ismap \\
            $options.args \\
            --t $task.cpus \\
            --output_dir $sampleName \\
            --queries $query_name \\
            --log ${prefix} \\
            --reference $reference_name \\
            --reads $reads

        # Reorganize output files
        mkdir ismapper
        mv $sampleName/*/* ismapper/


        # FIXME This would need to updated later for accommodating various static lengths per insertion element. Work with a CSV file name,element_sequence,sequence_length
        convert_ismapper_to_vcf.py

        """
}
