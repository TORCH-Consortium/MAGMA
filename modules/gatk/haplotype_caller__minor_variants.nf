process GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS {
    tag "$sampleName"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        tuple val(sampleName), path(bai), path(bam)
        path(reference)
        path("*")

    output:
        tuple val(sampleName), path("*.g.vcf.gz.tbi"), path("*.AllSites.g.vcf.gz"), emit: gvcf_ch
        tuple val(sampleName), path("*bam")


    script:

        """
        ${params.gatk_path} HaplotypeCaller --java-options "-Xmx${task.memory.giga}G" \\
            -R ${reference} \\
            -I ${bam} \\
            ${params.arguments} \\
            -bamout ${sampleName}.haplotype_caller.bam \\
            -O  ${sampleName}.AllSites.g.vcf.gz
        """

    stub:

        """
        echo "gatk HaplotypeCaller -Xmx${task.memory.giga}G \\
            -R ${reference} \\
            -I ${bam} \\
            ${params.arguments} \\
            -O  ${sampleName}.AllSites.g.vcf.gz"

        touch ${sampleName}.AllSites.g.vcf.gz
        """
}

