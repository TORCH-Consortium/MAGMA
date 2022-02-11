#!/usr/bin/env bash

set -xue

# NOTE: If there are problems in `mamba`, replace it with `conda``

mamba create -p ./xbs-nf-env-1 bioconda::gatk4=4.2.4.1 conda-forge::R=4.1  conda-forge::r-ggplot2=3.3.5 bioconda::datamash=1.1.0 bioconda::delly=0.8.7 bioconda::lofreq=2.1.5 bioconda::delly=0.8.7 bioconda::lofreq=2.1.5 bioconda::tb-profiler=4.1.1 bioconda::multiqc=1.11 bioconda::fastqc=0.11.8

mamba create -p ./xbs-nf-env-2 jemunro::quanttb=1.01 bioconda::bwa=0.7.17 bioconda::samtools=1.9 bioconda::iqtree=2.1.2 bioconda::snp-dists=0.8.2 bioconda::snp-sites=2.4.0 bioconda::bcftools=1.9 bioconda::snpeff=4.3.1t bioconda::clusterpicker=1.2.3
