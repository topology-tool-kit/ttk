#!/bin/bash

PARAVIEW_PATH=$1
PARAVIEW_VERSION_NEEDED=5.10.1

if [ -z "${PARAVIEW_PATH}" ];

then
    if ! [ -x "$(command -v paraview)" ]; then
        echo "No ParaView ${PARAVIEW_VERSION_NEEDED} installation found. Please provide path to ParaView ${PARAVIEW_VERSION_NEEDED} as first argument."
        echo "Usage:"
        echo "  $0 <Path to ParaView binary>"
        exit 0
    fi
    PV_VERSION_OUTPUT = $(paraview --version)
    PARAVIEW_VERSION = "${PV_VERSION_OUTPUT:17}"
    echo "$PARAVIEW_VERSION"
    
else
    PV_VERSION_OUTPUT=$($PARAVIEW_PATH --version)
    PARAVIEW_VERSION="${PV_VERSION_OUTPUT:17}"
    echo "$PARAVIEW_VERSION"
    if [[ "$PARAVIEW_VERSION" != "$PARAVIEW_VERSION_NEEDED" ]]; then
        echo "Path to ParaView $PARAVIEW_VERSION given, but ParaView $PARAVIEW_VERSION_NEEDED needed. Please provide path to correct version."
        exit 0
    fi

fi

DOCKER_ID=`docker run -d --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} ghcr.io/scivislab/ttk2:latest`
${PARAVIEW_PATH} --server-url=cs://127.0.0.1:11111 ${@:2}
docker kill ${DOCKER_ID} &> /dev/null
