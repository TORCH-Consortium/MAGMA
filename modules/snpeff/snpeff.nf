process SNPEFF {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(rawJointVariantsFile)
    path(ref_fasta)

    output:
    tuple val(joint_name), path("*.annotated.vcf")


    shell:

    '''
    gunzip -c !{rawJointVariantsFile} > !{joint_name}.temp.raw_variants.vcf

    sed -i 's/^!{ref_fasta.getBaseName()}/Chromosome/g' !{joint_name}.temp.raw_variants.vcf

    !{params.snpeff_path} \\
        !{params.arguments} \\
        !{joint_name}.temp.raw_variants.vcf \\
    > !{joint_name}.raw_variants.annotated.vcf

    rm !{joint_name}.temp.raw_variants.vcf

    sed -i 's/^Chromosome/!{ref_fasta.getBaseName()}/g' !{joint_name}.raw_variants.annotated.vcf
    '''

    stub:

    """
    touch ${joint_name}.raw_variants.annotated.vcf
    """

}
