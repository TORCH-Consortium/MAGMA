#!/bin/bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.

CONTAINER_TAG=0.4.0
DOCKER_NAMESPACE="ghcr.io/torch-consortium/magma"
CONTAINER_NAME="$DOCKER_NAMESPACE/biocontainer-ntmprofiler:$CONTAINER_TAG"
echo "Building container : $CONTAINER_NAME "

#============================================

#NTMDB_COMMIT="002fec58b78cd5bfbefd1106a1ac26777147ec6b"

#git clone https://github.com/ntm-db/ntm-db.git
#cd ntm-db
#git reset --hard  $NTMDB_COMMIT
#cd ../


#============================================

#NOTE: When changing the ntm-profiler version, update the 
#UPSTREAM_CONTAINER="quay.io/biocontainers/ntm-profiler:0.3.0--pyhdfd78af_0"
#docker pull $UPSTREAM_CONTAINER

docker build -t $CONTAINER_NAME .

#CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
#docker commit $CONTAINER_ID $CONTAINER_NAME
#docker push $DOCKER_NAMESPACE/$container_dir:$CONTAINER_TAG
#docker stop $CONTAINER_ID
