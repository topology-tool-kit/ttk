ZFP_VERSION=0.5.5

require-pkgs \
    vim             \
	curl		    \
	ca-certificates \
	cmake		    \
	build-essential	\
	ninja-build

fetch-src https://codeload.github.com/LLNL/zfp/tar.gz/0.5.5

cmake-default -DBUILD_EXAMPLES=OFF

