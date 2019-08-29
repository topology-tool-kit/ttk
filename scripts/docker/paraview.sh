#! /bin/bash

set -e

build_pkgs \
	build-essential		\
	curl			\
	cmake			\
	ninja-build		\
	python3.6-dev		\
	libexpat1-dev		\
	libeigen3-dev		\
	libfreetype6-dev	\
	liblz4-dev		\
	liblzma-dev		\
	libhdf5-dev		\
	libtiff-dev		\
	libjpeg-dev		\
	libpng-dev		\
	libogg-dev		\
	libtheora-dev		\
	libnetcdf-dev		\
	libnetcdf-cxx-legacy-dev\
	libxml2-dev		\
	libjsoncpp-dev		\
	libpugixml-dev		\
	libprotobuf-dev		\
	protobuf-compiler	\
	libcgns-dev		\
	libtbb-dev

runtime_pkgs \
	libstdc++6		\
	libpython3.6		\
	libexpat1		\
	libfreetype6		\
	liblz4-1		\
	liblzma5		\
	libhdf5-100		\
	libtiff5		\
	libjpeg8		\
	libpng16-16		\
	libogg0			\
	libtheora0		\
	libnetcdf13		\
	libnetcdf-c++4		\
	libxml++2.6-2v5		\
	libjsoncpp1		\
	libpugixml1v5		\
	libprotobuf10		\
	libcgns3.3		\
	libtbb2

BUILD_DIR=/root/paraview-build

echo "### build ParaView ###"

# get source code
mkdir -p $BUILD_DIR

curl -kL "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.6&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEW_VERSION}.tar.xz" | \
  tar Jx -C $BUILD_DIR --strip-components 1

# FIXME: hack to allow plugins built against ParaView to know about ospray targets 
echo "include(@ospray_DIR@/osprayConfig.cmake)" >> ${BUILD_DIR}/VTK/CMake/VTKConfig.cmake.in

# actually compile
mkdir -p $BUILD_DIR/build

pushd $BUILD_DIR/build

cmake -G Ninja \
    -DCMAKE_INSTALL_PREFIX=/usr \
    -DCMAKE_INSTALL_LIBDIR="lib" \
    -DCMAKE_BUILD_TYPE=Release \
    -DPARAVIEW_USE_VTKM=ON \
    -DPARAVIEW_BUILD_QT_GUI=OFF \
    -DPARAVIEW_ENABLE_PYTHON=ON \
    -DPARAVIEW_ENABLE_CATALYST=ON \
    -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
    -DVTK_USE_SYSTEM_LIBRARIES=ON \
    -DVTK_USE_SYSTEM_DOUBLECONVERSION=OFF \
    -DVTK_USE_SYSTEM_GL2PS:BOOL=OFF \
    -DVTK_USE_SYSTEM_PEGTL=OFF \
    -DVTK_USE_SYSTEM_XDMF2=OFF \
    -DVTK_USE_SYSTEM_GLEW=OFF \
    -DVTK_USE_X:BOOL=OFF \
    -DVTK_PYTHON_VERSION=3.6 \
    -DVTK_OPENGL_HAS_OSMESA:BOOL=ON \
    -DVTK_SMP_IMPLEMENTATION_TYPE=TBB \
    -DPARAVIEW_USE_OSPRAY=ON \
    -DOSPRAY_INSTALL_DIR=/usr \
    ..

ninja install

popd 

rm -rf $BUILD_DIR
