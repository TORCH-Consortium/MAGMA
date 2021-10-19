
process GATK_VARIANTS_TO_TABLE {
    tag "${joint_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    tuple val(joint_name), path(filteredSnpIncComplexVcfGz)

    output:
    path("*.95X.IncComplex.fa")


//FIXME
    shell:

    '''
    !{params.gatk_path} VariantsToTable -Xmx!{task.memory.giga}G \\
    -V !{filteredSnpIncComplexVcfGz} \\
    -GF GT \\
    -O /dev/stdout \\
    | sed -e 's/^\t//g' \\
    | sed -e 's/*/-/g' \\
    | sed -e 's/\./-/g' \\
    | sed '2,${/^.*\(-.*\)\{'"!{params.median_coverage_cutoff}"',\}.*$/d}' \\
    | !{params.datamash_path} transpose \\
    | sed -e 's/^/>/g' \\
    | sed -e 's/-GT/\n/g' \\
    | sed -e 's/\t//g' \\
    > !{joint_name}.95X.IncComplex.fa

    '''

    stub:

    """
    touch ${joint_name}.95X.IncComplex.fa
    """
}

