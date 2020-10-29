#! /bin/bash

set -e

build_pkgs \
    vim             \
	curl		    \
	ca-certificates \
	cmake		    \
	build-essential	\
	ninja-build

fetch_url https://codeload.github.com/LLNL/zfp/tar.gz/0.5.5 \
    | tar zx --strip-components 1

# actually compile
mkdir build 
pushd build

cmake -G Ninja \
      -DCMAKE_INSTALL_PREFIX=/usr	\
      -DCMAKE_BUILD_TYPE=Release	\
      -DBUILD_EXAMPLES=OFF 		\
      ..

cmake --build .
cmake --install .

popd

