#! /usr/bin/env bash
set -e

build_pkgs \
	build-essential	\
	pkg-config		\
	curl			\
	python3-dev		\
	python3-mako	\
	llvm-dev		\
    meson         	\
    ninja-build     \
	zlib1g-dev		\
	libdrm-dev		\
	gettext			\
	bison			\
	flex			\
	xz-utils

runtime_pkgs \
	libstdc++6		\
	libllvm10		\
	zlib1g

# get source

curl -kL https://mesa.freedesktop.org/archive/mesa-20.2.0.tar.xz \
    | tar Jx --strip-components 1

# configure and build
mkdir build

meson build 					\
    --prefix=/usr				\
    -Dosmesa=gallium			\
    -Dplatforms= 				\
    -Dgallium-drivers=swr,swrast\
    -Dglx=disabled				\
    -Dgles2=false				\
    -Dgles1=false				\
    -Dllvm=true					\
    -Ddri-drivers=				\
    -Dvulkan-drivers=			\
	-Dswr-arches=avx			\
    -Dshared-glapi=true

ninja -v -C build install

