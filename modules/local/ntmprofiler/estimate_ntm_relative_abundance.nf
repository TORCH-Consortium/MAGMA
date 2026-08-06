process NTMPROFILER_ESTIMATE_NTM_RELATIVE_ABUNDANCE {
    tag "${sampleName}"

    input:
    tuple val(sampleName), path(profileJson)

    output:
    tuple val(sampleName),
          path("${sampleName}.potential_NTM_fraction.txt"),
          emit: fraction

    tuple val(sampleName),
          path(profileJson),
          emit: enriched_json

    script:
    """
    ntmprofiler_json_extract_ntm_relative_abundance.py \
        ${profileJson} \
        ${sampleName}.potential_NTM_fraction.txt
    """

    stub:
    """
    echo "0" > ${sampleName}.potential_NTM_fraction.txt
    """
}
