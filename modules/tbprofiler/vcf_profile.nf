tb-profiler vcf_profile --db $RESISTANCE_DB -d $OUT_DIR/resistance/$JOINT_NAME/XBS $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP.RawIndels.vcf.gz


tb-profiler vcf_profile --lofreq_sample_name $i --db $RESISTANCE_DB -d $OUT_DIR/resistance/$JOINT_NAME/lofreq $OUT_DIR/minor_vars/$i.LoFreq.vcf.gz



process TBPROFILER_PROFILE {
    tag "${genomeName}"
    publishDir params.resultsDir_tbprofiler_per_sample, mode: params.saveMode_tbprofiler_per_sample, enabled: params.shouldPublish_tbprofiler_per_sample

    input:
    tuple val(genomeName), path(genomeReads)

    output:
    tuple path("results/*txt"), path("results/*json")


    script:
    """
    tb-profiler profile -1 ${genomeReads[0]} -2 ${genomeReads[1]}  -t ${task.cpus} -p $genomeName --txt
    """

    stub:
    """
    echo "tb-profiler profile -1 ${genomeReads[0]} -2 ${genomeReads[1]}  -t ${task.cpus} -p $genomeName --txt"

    mkdir results
    touch results/"${genomeName}.results.txt"
    touch results/"${genomeName}.results.json"
    """

}

