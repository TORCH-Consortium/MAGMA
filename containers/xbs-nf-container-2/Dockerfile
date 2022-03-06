FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com


RUN \
   micromamba install -y -n base  -c conda-forge -c bioconda \
     jemunro::quanttb=1.01 \
     bioconda::bwa=0.7.17 \
     bioconda::samtools=1.9 \
     bioconda::iqtree=2.1.2 \
     bioconda::snp-dists=0.8.2 \
     bioconda::snp-sites=2.4.0 \
     bioconda::bcftools=1.9 \
     bioconda::snpeff=4.3.1t \
     bioconda::clusterpicker=1.2.3 \
   && micromamba clean -a -y
