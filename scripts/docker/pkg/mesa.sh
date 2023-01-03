MESA_VERSION=20.3.5

# install packages required for build
require-pkgs \
    pkg-config      \
	python3-dev		\
	python3-mako	\
	llvm-dev		\
    meson         	\
	zlib1g-dev		\
	libdrm-dev		\
	gettext			\
	bison			\
	flex

# fetch and unpack source
fetch-src https://codeload.github.com/mesa3d/mesa/tar.gz/refs/tags/mesa-${MESA_VERSION}

# determine build type from CMake default
case ${CMAKE_BUILD_TYPE,,} in 
    relwithdebinfo)
        buildtype='debugoptimized'
        ;;
    minsizerel)
        buildtype='minsize'
        ;;
    *)
        buildtype=${CMAKE_BUILD_TYPE,,}
        ;;
esac

# configure 
meson build \
    -Dbuildtype=${buildtype}        \
    -Dosmesa=gallium			    \
    -Dplatforms= 				    \
    -Dgallium-drivers=swrast	    \
    -Dglx=disabled				    \
    -Dgles2=false				    \
    -Dgles1=false				    \
    -Dllvm=enabled				    \
    -Ddri-drivers=				    \
    -Dvulkan-drivers=			    \
    -Dswr-arches=				    \
    -Dshared-glapi=true

# build
ninja -C build 

# install
ninja -C build install

