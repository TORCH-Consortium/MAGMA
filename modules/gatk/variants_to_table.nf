nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/variants_to_table"
params.save_mode = 'copy'
params.should_publish = true


process GATK_VARIANTS_TO_TABLE {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    output:


    shell:

    '''
    gatk VariantsToTable -Xmx!{task.memory.giga}G \\
    -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP.ExDR.IncComplex.vcf.gz \\
    -GF GT \\
    -O /dev/stdout \\
    | sed -e 's/^\t//g' \\
    | sed -e 's/*/-/g' \\
    | sed -e 's/\./-/g' \\
    | sed '2,${/^.*\(-.*\)\{'"$CLUSTER_COVERAGE_THRESHOLD"',\}.*$/d}' \\
    | $DATAMASH transpose \\
    | sed -e 's/^/>/g' \\
    | sed -e 's/-GT/\n/g' \\
    | sed -e 's/\t//g' \\
    > $OUT_DIR/fasta/$JOINT_NAME/$JOINT_NAME.95X.IncComplex.fa

    '''

    stub:

    """

    """
}

