process TBPROFILER_VCF_PROFILE__COHORT {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(joint_name), path(mergedVcfIndex), path(mergedVcf)
        path(resistanceDb)

    output:
        path("results/*")


    script:
        def optionalDb  = resistanceDb ? "--db ${resistanceDb.name}" : ""

        """

        bcftools view ${mergedVcf} | sed 's/NC-000962-3-H37Rv/Chromosome/g' > intermediate.vcf

        cat  intermediate.vcf | bcftools view -Oz -o intermediate.vcf.gz

        ${params.tbprofiler_path} profile  \\
            ${optionalDb} \\
            --threads ${task.cpus}\\
            --vcf intermediate.vcf.gz \\
            ${params.arguments}


        """

    stub:
        """
        """

}
