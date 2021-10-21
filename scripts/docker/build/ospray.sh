#! /bin/bash

set -e

build_pkgs \
    build-essential         \
    curl                    \
    ca-certificates         \
    cmake                   \
    ninja-build             \
    libtbb-dev

runtime_pkgs \
    libtbb2

# --------------------------------------------------------------------------

case $PARAVIEW_VERSION in
5.[678]*)
    echo "using OSPRay 1.8.5"

    curl -L https://github.com/ospray/ospray/archive/v1.8.5.tar.gz | \
        tar xz --strip-components 1
    
    conf_args \
        -DOSPRAY_ENABLE_TUTORIALS=OFF           \
        -DOSPRAY_APPS_EXAMPLEVIEWER=OFF         \
        -DOSPRAY_APPS_UTILITIES=OFF             \
        -DOSPRAY_APPS_PARAVIEW_TFN_CVT=OFF      \
        -DOSPRAY_AUTO_DOWNLOAD_TEST_IMAGES=OFF  \
        -DOSPRAY_APPS_BENCHMARK=ON              \
        -DOSPRAY_ZIP_MODE=OFF
    ;;

*)
    echo "using OSPRay 2.x"

    curl -L https://github.com/ospray/ospray/archive/v2.4.0.tar.gz | \
        tar xz --strip-components 1

    conf_args \
        -DOSPRAY_ENABLE_APPS=OFF                \
        -DOSPRAY_MODULE_DENOISER=ON
    ;;
esac


mkdir build
pushd build

cmake -G Ninja \
    -DCMAKE_INSTALL_PREFIX=/usr             \
    -DCMAKE_BUILD_TYPE=Release	            \
    ${configure_args}                       \
    ..

cmake --build .
cmake --install .

popd

# REQUIRED - not sure why this isn't automatically called
ldconfig