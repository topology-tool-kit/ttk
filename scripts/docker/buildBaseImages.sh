#!/usr/bin/env bash

SCRIPTPATH=$(dirname $(realpath $0))

VERSIONS=(
    5.9.1-dev
    5.9.1-0.9.9

    5.9.0-dev
    5.9.0-0.9.9

    5.8.1-dev
    5.8.1-0.9.9

    5.8.0-dev
    5.8.0-0.9.9

    # no longer maintained
    # 5.7.0-stable
    # 5.6.2-0.9.8
    # 5.6.1-0.9.8
    # 5.6.0-0.9.8
)

for version in ${VERSIONS[@]}; do
    paraview=${version%-*}
    ttk=${version##*-}

    echo $paraview $ttk
    
    # build base image (without TTK)
    docker build \
        --build-arg paraview=${paraview} \
        -t topologytoolkit/ttk-base:${paraview} \
        -f build/Dockerfile.base \
        "${SCRIPTPATH}/build"

    # build dev image (without TTK)
    docker build \
       --build-arg dev=true \
       --build-arg paraview=${paraview} \
       -t topologytoolkit/ttk-dev:${paraview} \
       -f build/Dockerfile \
       "${SCRIPTPATH}/build"

    # build image (with TTK)
    docker build \
       --build-arg ttk=${ttk} \
       --build-arg paraview=${paraview} \
       -t topologytoolkit/ttk:${paraview}-${ttk} \
       -f build/Dockerfile \
       "${SCRIPTPATH}/build"
done

