FROM mambaorg/micromamba

# The conda env file has been copied via the build script
COPY --chown=$MAMBA_USER:$MAMBA_USER mtbseq-nf-env.yml /tmp/mtbseq-nf-env.yml

RUN micromamba install -y -f /tmp/mtbseq-nf-env.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  \
  && micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER GenomeAnalysisTK.jar /tmp/GenomeAnalysisTK.jar

RUN /opt/conda/bin/gatk-register /tmp/GenomeAnalysisTK.jar
