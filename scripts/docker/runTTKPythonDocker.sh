#!/bin/bash

docker run -it --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} ghcr.io/scivislab/ttk2:latest pvpython ${@:1}
