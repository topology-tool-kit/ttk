#! /bin/bash

set -e

echo "### build ZFP ###"

build_pkgs \
	build-essential	\
	curl		\
	ca-certificates \
	cmake		\
	ninja-build

BUILD_DIR=/root/zfp-build

# get source code
mkdir -p $BUILD_DIR

fetch_url https://codeload.github.com/LLNL/zfp/tar.gz/0.5.5 \
    | tar zx -C $BUILD_DIR --strip-components 1

# actually compile
mkdir -p $BUILD_DIR/build

pushd $BUILD_DIR/build

cmake -G Ninja \
      -DCMAKE_INSTALL_PREFIX=/usr	\
      -DCMAKE_BUILD_TYPE=Release 	\
      -DBUILD_EXAMPLES=OFF 		\
      ..

ninja install

popd

rm -rf $BUILD_DIR

