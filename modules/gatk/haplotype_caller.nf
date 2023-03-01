process GATK_HAPLOTYPE_CALLER {
    tag "$sampleName"
    label 'cpu_8_memory_8'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        tuple val(sampleName), path(bai), path(bam)
        path(ref_fasta)
        path("*")

    output:
        tuple val(sampleName), path("*.g.vcf.gz.tbi"), path("*.g.vcf.gz"), emit: gvcf_ch
        tuple val(sampleName), path("*bam")


    script:

        """
        ${params.gatk_path} HaplotypeCaller --java-options "-Xmx${task.memory.giga}G" \\
            -R ${ref_fasta} \\
            -I ${bam} \\
            -ERC GVCF \\
            ${params.arguments} \\
            -bamout ${sampleName}.haplotype_caller.bam \\
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

