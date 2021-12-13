process GATK_HAPLOTYPE_CALLER {
    tag "$sampleName"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        tuple val(sampleName), path(bai), path(bam)
        path(ref_fasta)
        path("*")

    output:
        tuple val(sampleName), path("*.g.vcf.gz.tbi"), path("*.g.vcf.gz")


    script:

        """
        ${params.gatk_path} HaplotypeCaller --java-options "-Xmx${task.memory.giga}G" \\
            -R ${ref_fasta} \\
            -I ${bam} \\
            -ERC GVCF \\
            ${params.arguments} \\
            -O ${sampleName}.g.vcf.gz
        """

    stub:

        """
        echo "gatk HaplotypeCaller -Xmx${task.memory.giga}G \\
            -R ${ref_fasta} \\
            -I ${bam} \\
            -ERC GVCF \\
            ${params.arguments} \\
            -O ${sampleName}.g.vcf.gz"

        touch ${sampleName}.g.vcf.gz
        touch ${sampleName}.g.vcf.gz.tbi
        """
}

