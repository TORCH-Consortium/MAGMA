
         $JAVA -Xmx64G -jar $GATK ApplyVQSR -R $REFERENCE -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_SNP.vcf.gz -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_inc-rRNA.vcf.gz --tranches-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.tranches --recal-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.recal.vcf.gz --ts-filter-level 99.90 -AS --exclude-filtered -mode SNP


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

