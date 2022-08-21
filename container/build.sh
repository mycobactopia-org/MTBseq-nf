#!/usr/bin/env bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.

DOCKER_NAMESPACE="rg.fr-par.scw.cloud/nfcontainers"

cp ../conda_envs/mtbseq-nf-env.yml ./
cp ../resources/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar ./

echo "Building mtbseq-nf container ..."
CONTAINER_TAG=0.9.5
CONTAINER_NAME=$DOCKER_NAMESPACE/"mtbseq-nf":$CONTAINER_TAG

echo "Container Name : $CONTAINER_NAME "
docker build -t $CONTAINER_NAME .
CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
docker commit $CONTAINER_ID $CONTAINER_NAME
docker push $DOCKER_NAMESPACE/"mtbseq-nf":$CONTAINER_TAG
docker stop $CONTAINER_ID


echo "Deleting the copied files"
rm mtbseq-nf-env.yml
rm GenomeAnalysisTK.jar
