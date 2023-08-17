FROM mambaorg/micromamba:0.25.1 AS base

#NOTE: The conda env file has been copied via the build script
COPY --chown=$MAMBA_USER:$MAMBA_USER magma-env-1.yml /tmp/magma-env-1.yml

RUN micromamba install -y -f /tmp/magma-env-1.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  && micromamba clean -a -y

FROM base AS condaenv

ENV PATH /opt/conda/bin:$PATH

COPY --chown=$MAMBA_USER:$MAMBA_USER resistance_db_who /tmp/who_catalogue_amended_v5.0.0

FROM condaenv AS resistanceDbWho

RUN cd /tmp/resistance_db_who && /opt/conda/bin/tb-profiler load_library resistance_db_who && cd /tmp

