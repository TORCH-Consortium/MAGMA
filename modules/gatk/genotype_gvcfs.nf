
$JAVA -Xmx64G -jar $GATK GenotypeGVCFs -R $REFERENCE -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.combined.vcf.gz -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation --sample-ploidy 1


nextflow.enable.dsl = 2


params.results_dir = "${params.outdir}/gatk4/haplotype_caller"
params.save_mode = 'copy'
params.should_publish = true


params.gatk_path = "gatk"
params.java_opts = "-Xmx4G"
params.contamination = 0

process GATK_HAPLOTYPE_CALLER {
    tag "${sampleId}_${interval_chunk_name}"
    label 'gatk4_container'

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    tuple val(sampleId),
            path(input_recal_merged_bam),
            path(input_recal_merged_bai),
            path(input_recal_merged_md5),
            val(scatter_id),
            val(interval_chunk_name),
            path(interval_list_file)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)


    output:

    tuple val(sampleId),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf"),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf.idx")


    script:

    """
    set -e

    ${params.gatk_path} --java-options "${params.java_opts}" \
                        HaplotypeCaller \
                        -R ${ref_fasta} \
                        -I ${input_recal_merged_bam} \
                        --output "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf" \
                        -contamination ${params.contamination} \
                        -ERC GVCF \
                        -L ${interval_list_file}
    """

    stub:

    """
    touch "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf" 
    touch "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf.idx" 

    """
}

