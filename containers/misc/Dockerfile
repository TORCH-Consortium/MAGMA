FROM mambaorg/micromamba:1.4.9 AS base

#NOTE: The conda env file has been copied via the build script
COPY --chown=$MAMBA_USER:$MAMBA_USER misc-env.yml /tmp/misc-env.yml

RUN micromamba install -y -f /tmp/misc-env.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  && micromamba clean -a -y

FROM base AS condaenv

ENV PATH /opt/conda/bin:$PATH
