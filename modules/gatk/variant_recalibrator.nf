
$JAVA -Xmx64G -jar $GATK VariantRecalibrator -R $REFERENCE -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_SNP.vcf.gz -AS --resource:coll2014,known=false,training=true,truth=true,prior=15.0 $RESOURCE_PATH/Coll2014.UVPapproved.rRNAexcluded.vcf.gz --resource:coll2018,known=false,training=true,truth=true,prior=15.0 $RESOURCE_PATH/Coll2018.UVPapproved.rRNAexcluded.vcf.gz --resource:Napier2020,known=false,training=true,truth=true,prior=15.0 $RESOURCE_PATH/Napier2020.UVPapproved.rRNAexcluded.vcf.gz --resource:Benavente2015,known=true,training=false,truth=false,prior=5.0 $RESOURCE_PATH/Benavente2015.UVPapproved.rRNAexcluded.vcf.gz -an AS_QD -an DP -an AS_MQ -mode SNP --output $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.recal.vcf.gz --tranches-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.tranches --target-titv 1.7 --truth-sensitivity-tranche 100.0 --truth-sensitivity-tranche 99.9 --truth-sensitivity-tranche 99.8 --truth-sensitivity-tranche 99.7 --truth-sensitivity-tranche 99.6 --truth-sensitivity-tranche 99.5 --truth-sensitivity-tranche 99.4 --truth-sensitivity-tranche 99.3 --truth-sensitivity-tranche 99.2 --truth-sensitivity-tranche 99.1 --truth-sensitivity-tranche 99.0 --max-gaussians 4 -mq-cap 60 --output-model $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.model --rscript-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.R #-an AS_MQRankSum #-an AS_FS -an AS_SOR # -an AS_ReadPosRankSum


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

