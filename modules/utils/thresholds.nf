// (XBS_merge#L14)

// this function checks whther the relative abundance threshold is met for one of the strains
function rel_abundance_threshold_met () {
    IFS=';'
    local i
    for i in ${@:2}
    do
        if [ 1 -eq "$(echo "$i >= $1" | bc)" ]
        then
            IFS=' '
            return 0
        fi
    done
    IFS=' '
    return 1
}


// this function checks whther the relative abundance threshold is met for one of the strains
function ntm_fraction_threshold_met () {
    IFS=';'
    local i
    for i in ${@:2}
    do
        if [ 1 -eq "$(echo "$i <= $1" | bc)" ]
        then
            IFS=' '
            return 0
        fi
    done
    IFS=' '
    return 1
}
