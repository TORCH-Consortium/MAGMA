/*
FIXME: Documentation comments

*/


process LOFREQ_CALL__NTM {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(recalibratedBam)
    path(reference)

    output:
    tuple val(sampleName), path("*.potential_NTM_fraction.txt")

    shell:

    '''
	lofreq call \\
	    -f !{reference} \\
	    -r !{reference.getBaseName()}:!{params.region} \\
        !{arguments} \\
	    !{recalibratedBam} \\
	| grep -v "#" \\
	| cut -f 2 -d ";" \\
	| tr -d 'AF=' \\
	| awk '{Total=Total+$1} END{print Total}' \\
	> !{sampleName}.potential_NTM_fraction.txt
    '''

    stub:

    """
    echo "${reference} -- ${reference.getBaseName()} -- ${params.region} -- ${sampleName} -- ${recalibratedBam}"

    echo "${params.arguments}"

    touch ${sampleName}.potential_NTM_fraction.txt
    """

}
