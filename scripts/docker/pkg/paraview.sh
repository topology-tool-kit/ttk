#! /bin/bash

set -e

require-pkgs \
    python3-dev		        \
    libexpat1-dev		    \
    libeigen3-dev		    \
    libfreetype6-dev	    \
    liblz4-dev			    \
    liblzma-dev			    \
    libhdf5-dev			    \
    libtiff-dev			    \
    libjpeg-dev			    \
    libpng-dev			    \
    libogg-dev			    \
    libtheora-dev		    \
    libnetcdf-dev		    \
    libnetcdf-cxx-legacy-dev\
    libxml2-dev			    \
    libjsoncpp-dev		    \
    libpugixml-dev		    \
    libprotobuf-dev		    \
    protobuf-compiler	    \
    libcgns-dev


RELEASE=${PARAVIEW_VERSION%.*}
fetch-src "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${RELEASE}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEW_VERSION}.tar.gz"

enable_rt=$(test -f /usr/local/lib/libospray.so && echo "ON" || echo "OFF")
enable_dn=$(test -f /usr/local/lib/libOpenImageDenoise.so && echo "ON" || echo "OFF")

export CMAKE_BUILD_TYPE=Release

conf_args \
    -DCMAKE_BUILD_TYPE=RELEASE                  \
    -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON     \
    -DPARAVIEW_ENABLE_RAYTRACING=${enable_rt}   \
    -DPARAVIEW_USE_VTKM=ON					    \
    -DPARAVIEW_USE_PYTHON=ON                    \
    -DPARAVIEW_USE_QT=OFF 				        \
    -DVTK_USE_X:BOOL=OFF 					    \
    -DVTK_OPENGL_HAS_OSMESA:BOOL=ON 		    \
    -DVTK_SMP_IMPLEMENTATION_TYPE=OpenMP        \
    -DVTKOSPRAY_ENABLE_DENOISER=${enable_dn}

cmake-default

# download and install materials
mkdir -p /usr/share/materials 
curl -kL https://gitlab.kitware.com/paraview/materials/-/archive/master/materials-master.tar.bz2 | \
    tar jxv --strip-components 1 -C /usr/share/materials
