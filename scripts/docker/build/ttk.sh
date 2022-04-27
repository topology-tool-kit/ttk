#! /bin/bash
set -e

build_pkgs \
    build-essential         \
    cmake                   \
    curl                    \
    libboost-system1.71-dev \
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
    libgraphviz-dev	        \
    libtbb-dev              \
    libtheora-dev           \
    libtiff-dev             \
    libxml2-dev             \
    ninja-build             \
    protobuf-compiler       \
    python3-dev             \
    python3-numpy-dev       \
    zlib1g-dev              \
    cmake-curses-gui        \
    git                     \
    gdb

runtime_pkgs \
    libboost-system1.71	\
    python3-numpy       \
    python3-sklearn     \
    libsqlite3-0        \
    graphviz            \
    libgomp1

if [ -n "${DEV}" ]; then
        echo "DEVELOPER MODE"
        exit
fi

# get source code
(curl -kL "https://github.com/topology-tool-kit/ttk/archive/${TTK_VERSION}.tar.gz" | tar zx --strip-components 1) ||
(curl -kL "https://github.com/topology-tool-kit/ttk/archive/v${TTK_VERSION}.tar.gz" | tar zx --strip-components 1)

# actually compile
mkdir build
pushd build

cmake -G Ninja \
    -DCMAKE_INSTALL_PREFIX=/usr \
    -DCMAKE_BUILD_TYPE=Release \
    -DTTK_BUILD_DOCUMENTATION=OFF \
    -DTTK_BUILD_PARAVIEW_PLUGINS=ON \
    -DTTK_BUILD_STANDALONE_APPS=OFF \
    -DTTK_BUILD_VTK_WRAPPERS=ON \
    -DTTK_BUILD_VTK_PYTHON_MODULE=OFF \
    -DTTK_ENABLE_DOUBLE_TEMPLATING=ON \
    -DTTK_ENABLE_CPU_OPTIMIZATION=OFF \
    -DTTK_ENABLE_OPENMP=ON \
    -DTTK_ENABLE_KAMIKAZE=ON \
    ..

# call Ninja manually to ignore duplicate targets
# cmake --build .
ninja -w dupbuild=warn install
cmake --install .

popd
