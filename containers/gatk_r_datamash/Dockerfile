FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com

RUN \
   micromamba install -y -n base -c defaults -c bioconda -c conda-forge \
      gatk4=4.2.0.0 \
      R=4.1 \
      r-ggplot2=3.3.5 \
      datamash=1.1.0 \
   && micromamba clean -a -y
