#!/bin/bash

PARAVIEW_VERSION=5.6.1

docker run -it --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} topologytoolkit/ttk:${PARAVIEW_VERSION}-master pvpython ${@:1}
