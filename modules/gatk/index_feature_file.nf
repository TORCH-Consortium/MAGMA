process GATK_INDEX_FEATURE_FILE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        tuple val(sampleName), path(vcf)
        val(outputPrefix)


    output:
        tuple val(sampleName), path("*.vcf.gz.tbi"), path(vcf), emit: sample_vcf_tuple
        tuple path("*.vcf.gz.tbi"), path(vcf), emit: vcf_tuple


    script:

        def outputFileArg = ( outputPrefix != "" ? "-O ${sampleName}.${outputPrefix}.vcf.gz.tbi" : "" )

        """
        ${params.gatk_path} IndexFeatureFile --java-options "-Xmx${task.memory.giga}G" \\
            ${outputFileArg} \\
            -I ${vcf}

        # Trigger the process again by changing this script
        """

    stub:

        def outputFileArg = ( outputPrefix != "" ? "-O ${sampleName}.${outputPrefix}.vcf.gz.tbi" : "" )

        """

        echo ${outputFileArg}

        touch ${sampleName}.${outputPrefix}.idx.vcf.gz
        touch ${sampleName}.${outputPrefix}.idx.vcf.gz.tbi

        """
}
