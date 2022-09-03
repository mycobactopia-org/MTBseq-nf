#!/usr/bin/env bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.
# NOTE: Login to github registry with `echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin`

CONTAINER_TAG=0.9.5

DOCKER_NAMESPACE="ghcr.io/mtb-bioinformatics"

echo "Downloading and uncompressing GATK jar"
wget "https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2"
tar -xf GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2 --wildcards '*.jar'
cp GenomeAnalysis*/*.jar .

echo "Coping the mtbseq-n env file"
cp ../conda_envs/mtbseq-nf-env.yml ./

echo "Building mtbseq-nf container ..."
CONTAINER_NAME=$DOCKER_NAMESPACE/"mtbseq-nf":$CONTAINER_TAG

echo "Container Name : $CONTAINER_NAME "
docker build -t $CONTAINER_NAME .
CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
docker commit "$CONTAINER_ID" "$CONTAINER_NAME"
docker push $DOCKER_NAMESPACE/"mtbseq-nf":$CONTAINER_TAG
docker stop "$CONTAINER_ID"


echo "Deleting the copied files"
rm mtbseq-nf-env.yml
rm -rf GenomeAnalysisTK*
