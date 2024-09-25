#!/bin/bash
##
## Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
##
## This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
##
## For quick overview of GPL-3 license, please refer
## https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
##
## - You MUST keep this license with original authors in your copy
## - You MUST acknowledge the original source of this software
## - You MUST state significant changes made to the original software
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program . If not, see <http://www.gnu.org/licenses/>.
##
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

CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
docker commit $CONTAINER_ID $CONTAINER_NAME
docker push $CONTAINER_NAME
docker stop $CONTAINER_ID
