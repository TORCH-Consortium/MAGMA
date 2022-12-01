FROM mambaorg/micromamba:0.25.1 AS base
MAINTAINER abhi18av@outlook.com

#NOTE: The conda env file has been copied via the build script
COPY --chown=$MAMBA_USER:$MAMBA_USER xbs-nf-env-1.yml /tmp/xbs-nf-env-1.yml

RUN micromamba install -y -f /tmp/xbs-nf-env-1.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  && micromamba clean -a -y

FROM base AS condaenv

ENV PATH /opt/conda/bin:$PATH

COPY --chown=$MAMBA_USER:$MAMBA_USER resistance_db_who /tmp/resistance_db_who

FROM condaenv AS resistanceDbWho

RUN cd /tmp/resistance_db_who && /opt/conda/bin/tb-profiler load_library resistance_db_who && cd /tmp

