process TBPROFILER_VCF_PROFILE__LOFREQ {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(name), path(mergedLofreqVcfIndex), path(mergedLofreqVcf)
        path(resistanceDb)

    output:
        tuple val(name), path("results/*")
        path("results/*"), emit: resistance_json


    script:
        def optionalDb  = resistanceDb ? "--db ${resistanceDb.name}" : ""

        """

        bcftools view ${mergedLofreqVcf} | sed 's/NC-000962-3-H37Rv/Chromosome/g' > intermediate.vcf
         
        cat  intermediate.vcf | bcftools view -Oz -o intermediate.vcf.gz

        ${params.tbprofiler_path} profile \\
            ${optionalDb} \\
            --threads ${task.cpus}\\
            --vcf intermediate.vcf.gz \\
            ${params.arguments}
        """

    stub:
        """
        mkdir results
        touch results/${sampleName}.results.json
        """

}
