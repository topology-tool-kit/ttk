#! /bin/bash

set -e

build_pkgs \
    build-essential		\
    curl				\
    ca-certificates		\
    cmake				\
    ninja-build			\
    python3.8-dev		\
    libexpat1-dev		\
    libeigen3-dev		\
    libfreetype6-dev	\
    liblz4-dev			\
    liblzma-dev			\
    libhdf5-dev			\
    libtiff-dev			\
    libjpeg-dev			\
    libpng-dev			\
    libogg-dev			\
    libtheora-dev		\
    libnetcdf-dev		\
    libnetcdf-cxx-legacy-dev\
    libxml2-dev			\
    libjsoncpp-dev		\
    libpugixml-dev		\
    libprotobuf-dev		\
    protobuf-compiler	\
    libcgns-dev			\
    libtbb-dev

runtime_pkgs \
    libstdc++6			\
    libpython3.8		\
    libexpat1			\
    libfreetype6		\
    liblz4-1			\
    liblzma5			\
    libhdf5-103			\
    libtiff5			\
    libjpeg8			\
    libpng16-16			\
    libogg0				\
    libtheora0			\
    libnetcdf15			\
    libnetcdf-c++4		\
    libxml++2.6-2v5		\
    libjsoncpp1			\
    libpugixml1v5		\
    libprotobuf17		\
    libcgns3.4			\
    libtbb2

RELEASE=${PARAVIEW_VERSION%.*}

# get source code
curl -kL "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${RELEASE}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEW_VERSION}.tar.xz" | \
  tar Jx --strip-components 1

# FIXME: hack to allow plugins built against ParaView to know about ospray targets
echo "include(@ospray_DIR@/osprayConfig.cmake)" >> VTK/CMake/VTKConfig.cmake.in

# actually compile
conf_args -G Ninja \
    -DCMAKE_INSTALL_PREFIX=/usr		 		\
    -DCMAKE_INSTALL_LIBDIR="lib" 			\
    -DCMAKE_BUILD_TYPE=Release 				\
    -DPARAVIEW_USE_VTKM=ON 					\
    -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
    -DVTK_USE_X:BOOL=OFF 					\
    -DVTK_PYTHON_VERSION=3.8 				\
    -DVTK_OPENGL_HAS_OSMESA:BOOL=ON 		\
    -DVTK_SMP_IMPLEMENTATION_TYPE=TBB

if [[ "$RELEASE" == "5.6" || "$RELEASE" == "5.7" ]] ; then
    echo "PARAVIEW 5.6/5.7 detected"

    sed -i -E 's/nullptr,(\s+\/\/ tp_)/0,      \1/g' VTK/Wrapping/PythonCore/PyVTK*.cxx
    sed -i -E 's/nullptr,(\s+\/\/ tp_)/0,      \1/g' VTK/Wrapping/Tools/vtkWrapPython*.c

    conf_args \
        -DPARAVIEW_BUILD_QT_GUI=OFF 		\
        -DPARAVIEW_ENABLE_PYTHON=ON 		\
        -DPARAVIEW_USE_OSPRAY=ON 			\
        -DOSPRAY_INSTALL_DIR=/usr 			\
        -DVTK_USE_SYSTEM_LIBRARIES=ON 		\
        -DVTK_USE_SYSTEM_DOUBLECONVERSION=OFF \
        -DVTK_USE_SYSTEM_GL2PS:BOOL=OFF 	\
        -DVTK_USE_SYSTEM_PEGTL=OFF 			\
        -DVTK_USE_SYSTEM_XDMF2=OFF 			\
        -DVTK_USE_SYSTEM_GLEW=OFF
else
    echo "PARAVIEW 5.8+ detected"
    conf_args \
        -DPARAVIEW_USE_QT=OFF 				\
        -DPARAVIEW_USE_PYTHON=ON 			\
        -DPARAVIEW_ENABLE_RAYTRACING=ON     \
        -DVTKOSPRAY_ENABLE_DENOISER=ON
fi

mkdir build
pushd build

cmake ${configure_args} ..
cmake --build .
cmake --install .

# download and install materials
mkdir -p /usr/share/materials 
curl -kL https://gitlab.kitware.com/paraview/materials/-/archive/master/materials-master.tar.bz2 | \
    tar jxv --strip-components 1 -C /usr/share/materials

popd
