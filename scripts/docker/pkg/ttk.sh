#! /bin/bash
set -e

require-pkgs \
    build-essential         \
    cmake                   \
    curl                    \
    libboost-system-dev     \
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
    libgraphviz-dev	    \
    libtheora-dev           \
    libtiff-dev             \
    libxml2-dev             \
    ninja-build             \
    protobuf-compiler       \
    python3-dev             \
    python3-numpy-dev       \
    zlib1g-dev
    
if [ -n "${DEV}" ]; then
        #echo "DEVELOPER MODE"
        exit
fi

# get source code
(curl -kL "https://github.com/topology-tool-kit/ttk/archive/${TTK_VERSION}.tar.gz" | tar zx --strip-components 1) ||
(curl -kL "https://github.com/topology-tool-kit/ttk/archive/v${TTK_VERSION}.tar.gz" | tar zx --strip-components 1)

# actually compile
cmake-default \
    -DTTK_BUILD_DOCUMENTATION=OFF \
    -DTTK_BUILD_PARAVIEW_PLUGINS=ON \
    -DTTK_BUILD_STANDALONE_APPS=OFF \
    -DTTK_BUILD_VTK_WRAPPERS=ON \
    -DTTK_BUILD_VTK_PYTHON_MODULE=OFF \
    -DTTK_ENABLE_DOUBLE_TEMPLATING=OFF \
    -DTTK_ENABLE_CPU_OPTIMIZATION=OFF \
    -DTTK_ENABLE_OPENMP=ON \
    -DTTK_ENABLE_KAMIKAZE=ON \
    ..

# call Ninja manually to ignore duplicate targets
# cmake --build .

# ninja -w dupbuild=warn install 
# cmake --install .

# popd
