FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com

RUN \
   micromamba install -y -n base -c defaults -c bioconda -c conda-forge \
    bwa=0.7.17 \
    samtools=1.9 \
   && micromamba clean -a -y
