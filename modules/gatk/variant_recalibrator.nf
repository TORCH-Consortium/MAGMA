nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/variant_recalibrator"
params.save_mode = 'copy'
params.should_publish = true


process GATK_VARIANT_RECALIBRATOR {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    output:

    script:

    """
    gatk VariantRecalibrator -Xmx${task.memory.giga}G \\
    -R $REFERENCE \\
    -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_SNP.vcf.gz \\
    --use-allele-specific-annotations \\
    --resource:coll2014,known=false,training=true,truth=true,prior=15.0 $RESOURCE_PATH/Coll2014.UVPapproved.rRNAexcluded.vcf.gz \\
    --resource:coll2018,known=false,training=true,truth=true,prior=15.0 $RESOURCE_PATH/Coll2018.UVPapproved.rRNAexcluded.vcf.gz \\
    --resource:Napier2020,known=false,training=true,truth=true,prior=15.0 $RESOURCE_PATH/Napier2020.UVPapproved.rRNAexcluded.vcf.gz \\
    --resource:Benavente2015,known=true,training=false,truth=false,prior=5.0 $RESOURCE_PATH/Benavente2015.UVPapproved.rRNAexcluded.vcf.gz \\
    -an AS_QD \\
    -an DP \\
    -an AS_MQ \\
    -mode SNP \\
    --tranches-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.tranches \\
    --target-titv 1.7 \\
    --truth-sensitivity-tranche 100.0 \\
    --truth-sensitivity-tranche 99.9 \\
    --truth-sensitivity-tranche 99.8 \\
    --truth-sensitivity-tranche 99.7 \\
    --truth-sensitivity-tranche 99.6 \\
    --truth-sensitivity-tranche 99.5 \\
    --truth-sensitivity-tranche 99.4 \\
    --truth-sensitivity-tranche 99.3 \\
    --truth-sensitivity-tranche 99.2 \\
    --truth-sensitivity-tranche 99.1 \\
    --truth-sensitivity-tranche 99.0 \\
    --max-gaussians 4 \\
    -mq-cap 60 \\
    --rscript-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.R \\
    #-an AS_MQRankSum \\
    #-an AS_FS \\
    # -an AS_SOR \\
    # -an AS_ReadPosRankSum \\
    --output $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.recal.vcf.gz \\
    --output-model $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.model

    """

    stub:

    """
    """
}
