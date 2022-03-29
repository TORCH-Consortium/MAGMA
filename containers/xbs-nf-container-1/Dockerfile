FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com

COPY --chown=$MAMBA_USER:$MAMBA_USER xbs-nf-env-1.yml /tmp/xbs-nf-env-1.yml

RUN micromamba install -y -f /tmp/xbs-nf-env-1.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  \
   && micromamba clean -a -y
