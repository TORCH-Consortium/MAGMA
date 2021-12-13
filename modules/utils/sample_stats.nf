process UTILS_SAMPLE_STATS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(samtoolsStats), path(wgsMetrics), path(flagStats), path(ntmFraction)

    output:
        path("*.stats.tsv")


    shell:

        '''

        # This function checks whether the NTM fraction threshold is met for one of the strains
        function ntm_fraction_threshold_met () {
            IFS=';'
            local i
            for i in ${@:2}
            do
                if [ 1 -eq "$(echo "$i <= $1" | bc)" ]
                then
                    IFS=' '
                    echo 1 # Return 1 because the sample passed the stats
                    return 0 # set a zero exitcode for if branching (i.e. this function passed)
                fi
            done
            IFS=' '
            echo 0 # Return 0 because the sample did not pass the stats
            return 1 # set a nonzero returncode to fail if statements
        }

        IFS=' '

        COVERAGE=$(cat !{wgsMetrics} | grep "^4411532" | cut -f 4)
        BREADTH_OF_COVERAGE=$(cat !{wgsMetrics} | grep "^4411532" | cut -f 14)
        TOTAL_SEQS=$(cat !{samtoolsStats} | grep "insert size average" | cut -f 3)
        MAPPED_P=$(cat !{flagStats} | grep "mapped (" | cut -f 2 -d "(" | cut -f 1 -d "%")
        INS_SIZE=$(cat !{samtoolsStats} | grep "raw total sequences" | cut -f 3)
        AVG_QUAL=$(cat !{samtoolsStats} | grep "average quality" | cut -f 3)
        WGS_METR=$(cat !{wgsMetrics} | grep "^4411532" | cut -f 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20,22,27)
        NTM_FRACTION=$(cat !{ntmFraction})

        if [ $COVERAGE -ge !{params.median_coverage_cutoff} ]
            then
                COVERAGE_THRESHOLD_MET=1
            else
                COVERAGE_THRESHOLD_MET=0
        fi


        if [ 1 -eq "$(echo "$BREADTH_OF_COVERAGE >= !{params.breadth_of_coverage_cutoff}" | bc)" ]
            then
                BREADTH_OF_COVERAGE_THRESHOLD_MET=1
            else
                BREADTH_OF_COVERAGE_THRESHOLD_MET=0
        fi


        if [ $COVERAGE -ge !{params.median_coverage_cutoff} ] && [ 1 -eq "$(echo "$BREADTH_OF_COVERAGE >= !{params.breadth_of_coverage_cutoff}" | bc)" ] && ntm_fraction_threshold_met !{params.ntm_fraction_cutoff} $NTM_FRACTION;
            then
                ALL_THRESHOLDS_MET=1
            else
                ALL_THRESHOLDS_MET=0
        fi


        echo -e "!{sampleName}\t${TOTAL_SEQS}\t${MAPPED_P}\t${INS_SIZE}\t${AVG_QUAL}\t${WGS_METR}\t${NTM_FRACTION}\t$(ntm_fraction_threshold_met !{params.ntm_fraction_cutoff} $NTM_FRACTION)\t${COVERAGE_THRESHOLD_MET}\t${BREADTH_OF_COVERAGE_THRESHOLD_MET}\t${ALL_THRESHOLDS_MET}" > !{sampleName}.stats.tsv
        '''
}
