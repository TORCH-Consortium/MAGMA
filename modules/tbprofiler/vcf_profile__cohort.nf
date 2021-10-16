

process TBPROFILER_VCF_PROFILE__SAMPLE {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(resistanceDb)

    output:


    script:
    """
    tb-profiler vcf_profile  \\
        --db $RESISTANCE_DB \\
        -d XBS \\
        $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP.RawIndels.vcf.gz

    """

    stub:
    """
    """

}
