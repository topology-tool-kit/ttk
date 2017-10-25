#!/bin/bash

# We need to build the ParaView or VTK dependencies from source
# on Travis CI since they're not provided through the container
# based infrastructure (I think an old VTK version is? might be too old)

MAKE_NUM_JOBS=2

if [ "$TTK_BUILD_PARAVIEW_PLUGINS" == "ON" ]; then
	wget "http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=source&os=all&downloadFile=ParaView-v5.4.0.tar.gz
	tar -xf ParaView-v5.4.0.tar.gz
	cd ttk/paraview/patch
	./patch-paraview-5.4.0.sh ../../../ParaView-v5.4.0/
	cd ../../../ParaView-v5.4.0/

	mkdir build
	cd build

	cmake -DCMAKE_INSTALL_PREFIX=../install \
		-DVTK_RENDERING_BACKEND=OpenGL2 \
		-DPARAVIEW_ENABLE_PYTHON=ON \
		-DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
		-DPARAVIEW_QT_VERSION=4 ..

	cmake --build . --target install -- -j $MAKE_NUM_JOBS
fi

