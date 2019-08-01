#!/bin/bash

PARAVIEW_PATH=$1
PARAVIEW_VERSION=5.6.1

if [ -z "${PARAVIEW_PATH}" ]; then
  echo "Usage:"
  echo "  $0 <Path to pvpython binary (version ${PARAVIEW_VERSION})>"
  exit 0
fi

docker run -it --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} topologytoolkit/ttk:${PARAVIEW_VERSION}-master pvpython ${@:2}
