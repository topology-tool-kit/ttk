#! /bin/bash

set -e

build_pkgs \
	build-essential \
	pkg-config	\
	curl		\
	python3-mako	\
	python3-dev	\
	python3-mako	\
	llvm-dev	\
        meson           \
        ninja-build     \
	zlib1g-dev	\
	libdrm-dev	\
	gettext		\
	bison		\
	flex		\
	xz-utils

runtime_pkgs \
	libstdc++6	\
	libllvm6.0	\
	zlib1g		

BUILD_DIR=/root/mesa-build

# get source
mkdir -p $BUILD_DIR

curl -kL https://mesa.freedesktop.org/archive/mesa-19.1.0.tar.xz \
    | tar Jx -C $BUILD_DIR --strip-components 1

# configure and build
pushd $BUILD_DIR

meson build/ 				\
    --prefix=/usr			\
    -Dosmesa=gallium			\
    -Dplatforms= 			\
    -Dgallium-drivers=swrast,swr	\
    -Dglx=disabled			\
    -Dgles2=false			\
    -Dgles1=false			\
    -Dllvm=true				\
    -Ddri-drivers=			\
    -Dvulkan-drivers=			\
    -Dshared-glapi=true

ninja -C build install

popd

rm -rf $BUILD_DIR
