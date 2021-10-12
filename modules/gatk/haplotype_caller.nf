nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/haplotype_caller"
params.save_mode = 'copy'
params.should_publish = true
params.arguments = " -ploidy 1 --read-filter MappingQualityNotZeroReadFilter -G StandardAnnotation -G AS_StandardAnnotation"


process GATK_HAPLOTYPE_CALLER {
    tag "$sampleName"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    tuple val(sampleName), path(bai), path(bam)
    path(ref_fasta)

    output:
    tuple val(sampleName), path("*.g.vcf.gz")


    script:

    """
    gatk HaplotypeCaller -Xmx${task.memory.giga}G \\
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
    """
}

