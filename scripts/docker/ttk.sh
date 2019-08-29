#! /bin/bash

set -e

BUILD_DIR=/root/ttk-build

build_pkgs \
        build-essential         \
        cmake                   \
        curl                    \
        libboost-system1.62-dev \
        libcgns-dev             \
        libeigen3-dev           \
        libexpat1-dev           \
        libfreetype6-dev        \
        libhdf5-dev             \
        libjpeg-dev             \
        libjsoncpp-dev          \
        liblz4-dev              \
        liblzma-dev             \
        libnetcdf-cxx-legacy-dev\
        libnetcdf-dev           \
        libogg-dev              \
        libpng-dev              \
        libprotobuf-dev         \
        libpugixml-dev          \
        libsqlite3-dev          \
	libgraphviz-dev		\
        libtbb-dev              \
        libtheora-dev           \
        libtiff-dev             \
        libxml2-dev             \
        ninja-build             \
        protobuf-compiler       \
        python3-dev             \
        python3-numpy-dev       \
        zlib1g-dev

runtime_pkgs \
	libboost-system1.62	\
	python3-numpy		\
	python3-sklearn		\
	libsqlite3-0		\
	graphviz		\
	libgomp1

echo "### build TTK ###"

# get source code
mkdir -p $BUILD_DIR

curl -kL "https://github.com/topology-tool-kit/ttk/archive/${TTK_VERSION}.tar.gz" | \
  tar zx -C $BUILD_DIR --strip-components 1

# actually compile
mkdir -p $BUILD_DIR/build

pushd $BUILD_DIR/build

cmake -G Ninja \
    -DCMAKE_INSTALL_PREFIX=/usr \
    -DCMAKE_BUILD_TYPE=Release \
    -DTTK_BUILD_DOCUMENTATION=OFF \
    -DTTK_BUILD_PARAVIEW_PLUGINS=ON \
    -DTTK_BUILD_STANDALONE_APPS=OFF \
    -DTTK_BUILD_VTK_WRAPPERS=ON \
    -DVTK_DIR=/usr/lib/cmake/paraview-5.6 \
    -DTTK_INSTALL_PLUGIN_DIR=${PV_PLUGIN_PATH} \
    ..

ninja install

popd 

rm -rf $BUILD_DIR
