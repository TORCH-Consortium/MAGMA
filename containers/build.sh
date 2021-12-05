#!/usr/bin/env bash

#!/bin/bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.

DOCKER_NAMESPACE="abhi18av"

for container_dir in $(find * -type d); do
  echo "Building $container_dir ..."
  cd $container_dir
  CONTAINER_TAG=0.0.1
  CONTAINER_NAME=$DOCKER_NAMESPACE/$container_dir:$CONTAINER_TAG
  docker build -t $CONTAINER_NAME .
  CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
  docker commit $CONTAINER_ID $CONTAINER_NAME
  docker push quay.io/$DOCKER_NAMESPACE/$container_dir
  docker stop $CONTAINER_ID
  cd ..
done
