#!/bin/bash

# We need to build the ParaView or VTK dependencies from source
# on Travis CI since they're not provided through the container
# based infrastructure (I think an old VTK version is? might be too old)

if [ "$TTK_BUILD_PARAVIEW_PLUGINS" == "ON" ]; then
	wget "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=binary&os=Sources&downloadFile=ParaView-v5.4.1.tar.gz" -O ParaView-v5.4.1.tar.gz
	tar -xf ParaView-v5.4.1.tar.gz
	cd ParaView-v5.4.1

	mkdir build
	cd build

	cmake -DCMAKE_INSTALL_PREFIX=../install \
		-DVTK_RENDERING_BACKEND=OpenGL2 \
		-DPARAVIEW_ENABLE_PYTHON=ON \
		-DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
		-DPARAVIEW_QT_VERSION=4 ..

	cmake --build . --target install -- -j `nproc`

elif [ "$TTK_BUILD_VTK_WRAPPERS" == "ON" ]; then
	wget "http://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz"
	tar -xf VTK-7.1.1.tar.gz
	cd VTK-7.1.1

	cmake -DCMAKE_INSTALL_PREFIX=./install \
		-DVTK_RENDERING_BACKEND=OpenGL2 \
		-DVTK_WRAP_PYTHON=ON .

	cmake --build . --target install -- -j `nproc`

fi

