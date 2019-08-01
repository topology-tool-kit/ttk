#!/bin/bash

PARAVIEW_VERSION=5.6.1

if [ -z "${PARAVIEW_PATH}" ]; then
  echo "Usage:"
  echo "  $0 [<Standard pvpython arguments (Python script, data, etc.)>]"
  exit 0
fi

docker run -it --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} topologytoolkit/ttk:${PARAVIEW_VERSION}-master pvpython ${@:1}
