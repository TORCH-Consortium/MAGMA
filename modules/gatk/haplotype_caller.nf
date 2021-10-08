nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/haplotype_caller"
params.save_mode = 'copy'
params.should_publish = true


process GATK_HAPLOTYPE_CALLER {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:


    output:


    script:

    """
    gatk HaplotypeCaller -Xmx${task.memory.giga}G \\
	-R $REFERENCE \\
    -I $OUT_DIR/mapped/$SAMPLE_ID.recal_reads.bam \\
    -ploidy 1 \\
    -ERC GVCF \\
    --read-filter MappingQualityNotZeroReadFilter \\
    -G StandardAnnotation \\
    -G AS_StandardAnnotation \\
    -O $OUT_DIR/gvcf/$SAMPLE_ID.g.vcf.gz
    """

    stub:

    """

    """
}

