EMBREE_VERSION=3.13.2

require-pkgs \
    libtbb-dev

# get source code
fetch-src https://github.com/embree/embree/archive/v${EMBREE_VERSION}.tar.gz

cmake-default \
    -DEMBREE_TASKING_SYSTEM=TBB \
    -DEMBREE_TUTORIALS=OFF

