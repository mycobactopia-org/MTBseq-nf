FROM mambaorg/micromamba:0.25.1 AS base

#NOTE: The conda env file has been copied via the build script
COPY --chown=$MAMBA_USER:$MAMBA_USER mtbseq-nf-env.yml /tmp/mtbseq-nf-env.yml

RUN micromamba create -n mtbseq-nf-env

RUN micromamba install -y -n mtbseq-nf-env -f /tmp/mtbseq-nf-env.yml && micromamba clean -a -y

FROM base AS condaenv

ENV PATH /opt/conda/envs/mtbseq-nf-env/bin:$PATH

COPY --chown=$MAMBA_USER:$MAMBA_USER GenomeAnalysisTK.jar /tmp/GenomeAnalysisTK.jar

FROM condaenv AS gatkjar

RUN /opt/conda/envs/mtbseq-nf-env/bin/gatk-register /tmp/GenomeAnalysisTK.jar
