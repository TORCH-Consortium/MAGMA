/*
FIXME: Documentation comments

*/



process GATK_VARIANT_RECALIBRATOR {
    tag "${joint_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    output:


    //TODO: A good candidate to refactor later on to accomodate dynamic number of files

    script:

    """
    gatk VariantRecalibrator -Xmx${task.memory.giga}G \\
        -R ${reference} \\
        -V ${raw_snp} \\
        --resource:coll2014,known=false,training=true,truth=true,prior=15.0 ${coll2014_vcf} \\
        --resource:coll2018,known=false,training=true,truth=true,prior=15.0 ${coll2018_vcf} \\
        --resource:Napier2020,known=false,training=true,truth=true,prior=15.0 ${napier2020_vcf} \\
        --resource:Benavente2015,known=true,training=false,truth=false,prior=5.0 ${benavente2015} \\
        ${params.arguments} \\
        -mode ${params.mode} \\
        --tranches-file ${joint_name}.snp.tranches \\
        --rscript-file ${joint_name}.snp.R \\
        --output ${joint_name}.snp.recal.vcf.gz \\
        --output-model ${joint_name}.snp.model
    """

    stub:

    """
    """
}
