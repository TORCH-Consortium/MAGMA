#!/bin/bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.

CONTAINER_TAG=1.2.0
DOCKER_NAMESPACE="ghcr.io/torch-consortium/magma"

CONTAINER_NAME="$DOCKER_NAMESPACE/magma-container-1:$CONTAINER_TAG"

echo "Building container : $CONTAINER_NAME "

cp ../../conda_envs/magma-env-1.yml ./


docker build -t $CONTAINER_NAME .
CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
docker commit $CONTAINER_ID $CONTAINER_NAME
docker push $DOCKER_NAMESPACE/$container_dir:$CONTAINER_TAG
docker stop $CONTAINER_ID
