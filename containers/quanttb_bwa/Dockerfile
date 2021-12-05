FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com

RUN \
   micromamba install -y -n base -c defaults -c bioconda -c conda-forge -c jemunro \
    quanttb=1.01 \
    bwa=0.7.17 \
   && micromamba clean -a -y
