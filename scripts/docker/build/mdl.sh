#! /bin/bash

set -e

build_pkgs \
	build-essential 		\
	pkg-config				\
	curl					\
	cmake					\
	cmake-curses-gui		\
    ninja-build     		\
	libfreeimage-dev		\
	python3-minimal			\
	clang-7					\

runtime_pkgs \
	libstdc++6				\
	libfreeimage3			\
	zlib1g


# get source code
curl -kL https://github.com/NVIDIA/MDL-SDK/archive/2020.1.2.tar.gz | \
	tar xz --strip-components 1

mkdir build
pushd build

exit 

cmake -G Ninja \
    -DCMAKE_BUILD_TYPE=Release    \
    -DCMAKE_INSTALL_PREFIX=/usr   \
    -DBUILD_TESTING=OFF           \
    ..

cmake --build .
cmake --install .

popd

exit 0

