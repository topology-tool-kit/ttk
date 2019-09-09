#!/bin/bash

PARAVIEW_PATH=$1
PARAVIEW_VERSION=5.6.1

if [ -z "${PARAVIEW_PATH}" ]; then
  echo "Usage:"
  echo "  $0 <Path to ParaView binary (version ${PARAVIEW_VERSION})>"
  exit 0
fi

DOCKER_ID=`docker run -d --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} topologytoolkit/ttk:${PARAVIEW_VERSION}-master`
${PARAVIEW_PATH} --server-url=cs://127.0.0.1:11111 ${@:2}
docker kill ${DOCKER_ID} &> /dev/null
