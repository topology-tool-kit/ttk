#! /bin/bash

set -e

case $PARAVIEW_VERSION in
5.[678]*)
    ;;
    
*)
    build_pkgs \
        build-essential         \
        curl                    \
        ca-certificates         \
        cmake                   \
        ninja-build             \
        python3-minimal         \
        libtbb-dev

    runtime_pkgs \
        libtbb2

    # --------------------------------------------------------------------------

    # get source code
    curl -L https://github.com/OpenImageDenoise/oidn/releases/download/v1.2.4/oidn-1.2.4.src.tar.gz | \
    tar xz --strip-components 1

    mkdir build
    pushd build

    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release    \
        -DCMAKE_INSTALL_PREFIX=/usr   \
        -DOIDN_APPS=OFF               \
        ..

    cmake --build .
    cmake --install .

    popd
    ;;
esac