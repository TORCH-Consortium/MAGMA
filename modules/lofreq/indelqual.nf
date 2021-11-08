process LOFREQ_INDELQUAL {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(recalibratedBam)
    path(ref_fasta)

    output:
    tuple val(sampleName), path("*.dindel.bam")

    script:

    """
    ${params.lofreq_path} indelqual \\
        -f ${ref_fasta} \\
        --dindel \\
        -o ${sampleName}.dindel.bam \\
        ${recalibratedBam}
    """

    stub:

    """
    echo "lofreq indelqual \\
        -f ${ref_fasta} \\
        --dindel \\
        -o ${sampleName}.dindel.bam \\
        ${recalibratedBam} "

    touch ${sampleName}.potential_NTM_fraction.txt
    touch ${sampleName}.dindel.bam
    """

}
