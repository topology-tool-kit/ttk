#!/bin/bash

# We need to build the ParaView or VTK dependencies from source
# on Travis CI since they're not provided through the container
# based infrastructure (I think an old VTK version is? might be too old)

MAKE_NUM_JOBS=2

# assert pwd == <ttk-folder>/ext

# TODO: may want to install VTK only
# if TTK_BUILD_PARAVIEW_PLUGINS is OFF
if [ "$TTK_BUILD_VTK_WRAPPERS" == "ON" ]; then
  mkdir ParaView-nightly
  cd ParaView-nightly || exit 2
  git init
  git remote add origin "https://gitlab.kitware.com/paraview/paraview.git"
  git fetch --depth 1 origin "864d18d1b5500cb4088ecdc56d2ffd9922c81121"
  git checkout FETCH_HEAD
  git submodule update --init --recursive
  cd ..

  cd ../paraview/patch || exit 2
  ./patch-paraview-5.7.0.sh ../../ext/ParaView-nightly/
  cd ../../ext/ParaView-nightly/ || exit 2

  mkdir build
  cd build || exit 2

  cmake -DCMAKE_INSTALL_PREFIX=../install \
    -DVTK_RENDERING_BACKEND=OpenGL2 \
    -DPARAVIEW_ENABLE_PYTHON=ON \
    -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
    -DPARAVIEW_QT_VERSION=4 ..

  cmake --build . --target install -- -j $MAKE_NUM_JOBS
fi

