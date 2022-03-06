#!/bin/bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.

DOCKER_NAMESPACE="rg.nl-ams.scw.cloud/xbs-nf-containers"

for container_dir in $(find * -type d); do
  echo "Building $container_dir ..."
  cd $container_dir
  CONTAINER_TAG=0.1.0
  CONTAINER_NAME=$DOCKER_NAMESPACE/$container_dir:$CONTAINER_TAG
  echo "Container Name : $CONTAINER_NAME "
  docker build -t $CONTAINER_NAME .
  CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
  docker commit $CONTAINER_ID $CONTAINER_NAME
  docker push $DOCKER_NAMESPACE/$container_dir:$CONTAINER_TAG
  docker stop $CONTAINER_ID
  cd ..
done
