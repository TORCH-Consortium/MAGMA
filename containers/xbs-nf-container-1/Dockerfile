FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com


RUN \
   micromamba install -y -n base  -c conda-forge -c bioconda \
	bioconda::gatk4=4.2.4.1 \
	conda-forge::R=4.1 \
	conda-forge::r-ggplot2=3.3.5 \
	bioconda::datamash=1.1.0 \
	bioconda::delly=0.8.7 \
	bioconda::lofreq=2.1.5 \
	bioconda::delly=0.8.7 \
	bioconda::lofreq=2.1.5 \
	bioconda::tb-profiler=4.1.1 \
	bioconda::multiqc=1.11 \
	bioconda::fastqc=0.11.8 \
   && micromamba clean -a -y


RUN micromamba install -y -n base  conda-forge::procps-ng  \
   && micromamba clean -a -y
