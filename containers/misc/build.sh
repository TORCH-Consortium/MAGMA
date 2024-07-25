#!/bin/bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.

CONTAINER_TAG=2.0.0-theta
DOCKER_NAMESPACE="ghcr.io/torch-consortium/magma"
CONTAINER_DIR=misc

CONTAINER_NAME=$DOCKER_NAMESPACE/$CONTAINER_DIR:$CONTAINER_TAG

echo "Building container : $CONTAINER_NAME "


docker build -t $CONTAINER_NAME .
CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
docker commit $CONTAINER_ID $CONTAINER_NAME
docker push $DOCKER_NAMESPACE/$CONTAINER_DIR:$CONTAINER_TAG
docker stop $CONTAINER_ID
