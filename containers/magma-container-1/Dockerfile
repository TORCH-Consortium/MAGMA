FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com

COPY --chown=$MAMBA_USER:$MAMBA_USER magma-env-1.yml /tmp/magma-env-1.yml

RUN micromamba install -y -f /tmp/magma-env-1.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  \
   && micromamba clean -a -y
