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

# install ISPC
curl -L https://github.com/ispc/ispc/releases/download/v1.12.0/ispc-v1.12.0b-linux.tar.gz | \
  tar xz -C /usr/bin --strip-components 2 --wildcards "*/bin/ispc"

# --------------------------------------------------------------------------

# get source code
curl -L https://github.com/embree/embree/archive/v3.8.0.tar.gz | \
  tar xz --strip-components 1

mkdir build
pushd build

cmake -G Ninja \
  -DCMAKE_BUILD_TYPE=Release	\
  -DEMBREE_TUTORIALS=OFF        \
  -DEMBREE_ISA_SSE2=ON          \
  -DEMBREE_ISA_SSE42=ON         \
  -DEMBREE_ISA_AVX=ON           \
  -DEMBREE_ISA_AVX2=ON          \
  -DEMBREE_ISA_AVX512SKX=OFF    \
  -DEMBREE_ISA_AVX512KNL=OFF    \
  ..

cmake --build .
cmake --install .

popd
