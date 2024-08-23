process SNPEFF {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(joint_name), path(rawJointVariantsFile)
        path(ref_fasta)

    output:
        tuple val(joint_name), path("*.snpeff.vcf")


    shell:

        '''
        rename_vcf_chrom.py --vcf !{rawJointVariantsFile}  --source !{params.ref_fasta_basename} --target 'Chromosome' \
            | !{params.snpeff_path} -nostats !{params.arguments}  \
            | rename_vcf_chrom.py --target !{params.ref_fasta_basename} --source 'Chromosome' \
         > !{joint_name}.temp.vcf

         cp {joint_name}.temp.vcf  {joint_name}.raw_variants.snpeff.vcf
        '''

    stub:

        """
        touch ${joint_name}.raw_variants.annotated.vcf
        """

}
