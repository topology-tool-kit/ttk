#!/usr/bin/env bash

SCRIPTPATH=$(dirname $(realpath $0))

VERSIONS=(
    5.8.1-dev
    5.8.0-dev
    5.7.0-dev
    5.6.2-stable
    5.6.1-stable
    5.6.0-stable
)

for version in ${VERSIONS[@]}; do
    paraview=${version%-*}
    ttk=${version#*-}

    # build base image (without TTK)
    docker build \
        --build-arg paraview=${paraview} \
        -t topologytoolkit/ttk-base:${paraview} \
        -f build/Dockerfile.base \
        "${SCRIPTPATH}/build"

    # build dev image (with TTK)
    docker build \
        --build-arg ttk=${ttk} \
        --build-arg paraview=${paraview} \
        -t topologytoolkit/ttk:${paraview}-${ttk} \
        -f build/Dockerfile \
        "${SCRIPTPATH}/build"
done

# docker build --build-arg ttk=dev --build-arg paraview=5.8.1 -t topologytoolkit/ttk:5.8.1-dev -f build/Dockerfile build

