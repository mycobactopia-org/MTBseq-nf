FROM mambaorg/micromamba
MAINTAINER Abhinav Sharma <abhi18av@outlook.com>

COPY --chown=$MAMBA_USER:$MAMBA_USER mtbseq-nf-env.yml /tmp/mtbseq-nf-env.yml

RUN micromamba install -y -f /tmp/mtbseq-nf-env.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  \
   && micromamba clean -a -y
