#!/usr/bin/env bash

set -e

# NOTE: Please replace `conda` with `mamba` if it is installed for faster installs.

condaBinary="conda"

#NOTE: The conda environments are expected by the `conda_local` profile to be created within `conda_envs` directory

#NOTE: Adding a step to automatically register the gatk jar
echo "Downloading GATK jar"

wget "https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2"

tar -xf GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2 --wildcards '*.jar'

cp GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar .

echo "Creating mtbseq-nf-env"
$condaBinary env create -p mtbseq-nf-env --file mtbseq-nf-env.yml

echo "Registering GATK Jar"

# TODO: Good candidate for a clean approach in a refactor.
eval "$(conda shell.bash hook)"
$condaBinary activate "./mtbseq-nf-env"

gatk-register GenomeAnalysisTK.jar 

echo "Testing mtbseq-nf-env"
MTBseq --version
MTBseq --check

echo "Cleaning files"
rm -rf GenomeAnalysisTK*

echo "Cleaning files"
rm -rf GenomeAnalysisTK*
