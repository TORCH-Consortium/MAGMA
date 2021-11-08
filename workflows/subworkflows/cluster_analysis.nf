

workflow CLUSTER_ANALYSIS {

    take:
        incComplexFiles //tuple path(incComplexFasta), path(incComplexNewickTree)
        exComplexFile //tuple path(exComplexFasta), path(exComplexNewickTree)


    main:
        CLUSTER_PICKER(incComplexFiles, 5)
        CLUSTER_PICKER(incComplexFiles, 12)
        CLUSTER_PICKER(exComplexFiles, 5)
        CLUSTER_PICKER(exComplexFiles, 12)
}
