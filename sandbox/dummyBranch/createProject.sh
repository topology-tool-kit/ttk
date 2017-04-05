#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

# ex: ScalarFieldSmoother
Name=$1 

smallName="`tr '[:upper:]' '[:lower:]'  <<< ${Name:0:1}`${Name:1}"
bigName=`echo $Name | tr '[:lower:]' '[:upper:]'`

if [ -z "${Name}" ]; then
  echo "Usage:"
  echo "  $0 <Name, first letter in uppercase, no space>"
  echo "Example:"
  echo "  $0 ScalarFieldSmoother"
  exit -1
fi

# check for dependencies
SED=`which sed 2> /dev/null`
if [ -z "$SED" ]; then
  echo "Error: Please install sed."
  exit -1
fi


echo "Creating project ${smallName} ${Name} ${bigName}..."

function replace {
  $SED "s/blank/${smallName}/g" $1 > tmp.cmake
  mv tmp.cmake $1
  $SED "s/Blank/${Name}/g" $1 > tmp.cmake
  mv tmp.cmake $1
  $SED "s/BLANK/${bigName}/g" $1 > tmp.cmake
  mv tmp.cmake $1
}

mkdir -p sandbox/baseCode 2> /dev/null
mkdir -p sandbox/vtkWrappers 2> /dev/null
mkdir -p standalone 2> /dev/null
mkdir -p paraview/client 2> /dev/null
mkdir -p paraview/server 2> /dev/null

# 1) duplicate the blank base code
cp -R ../../core/baseCode/blank sandbox/baseCode/${smallName}
mv sandbox/baseCode/${smallName}/blank.cmake \
  sandbox/baseCode/${smallName}/${smallName}.cmake
replace "sandbox/baseCode/${smallName}/${smallName}.cmake"
mv sandbox/baseCode/${smallName}/Blank.cpp \
  sandbox/baseCode/${smallName}/${Name}.cpp
replace "sandbox/baseCode/${smallName}/${Name}.cpp"
mv sandbox/baseCode/${smallName}/Blank.h \
  sandbox/baseCode/${smallName}/${Name}.h
replace "sandbox/baseCode/${smallName}/${Name}.h"

# 2) duplicate the blank wrapper
cp -R ../../core/vtkWrappers/vtkBlank sandbox/vtkWrappers/vtk${Name}
mv sandbox/vtkWrappers/vtk${Name}/vtkBlank.cmake \
  sandbox/vtkWrappers/vtk${Name}/vtk${Name}.cmake
replace "sandbox/vtkWrappers/vtk${Name}/vtk${Name}.cmake"
mv sandbox/vtkWrappers/vtk${Name}/vtkBlank.cpp \
  sandbox/vtkWrappers/vtk${Name}/vtk${Name}.cpp
replace "sandbox/vtkWrappers/vtk${Name}/vtk${Name}.cpp"
mv sandbox/vtkWrappers/vtk${Name}/vtkBlank.h \
  sandbox/vtkWrappers/vtk${Name}/vtk${Name}.h
replace "sandbox/vtkWrappers/vtk${Name}/vtk${Name}.h"

# 3) duplicate the blank standalone projects
mkdir standalone/${smallName}
cp -R ../../standalone/Blank/cmd/ standalone/${smallName}/cmd
rm -R standalone/${smallName}/cmd/build &> /dev/null
replace "standalone/${smallName}/cmd/CMakeLists.txt"
replace "standalone/${smallName}/cmd/main.cpp"
rm standalone/${smallName}/cmd/core &> /dev/null
ln -sf ../../../../../core standalone/${smallName}/cmd/core
ln -sf ../../../sandbox standalone/${smallName}/cmd/sandbox

cp -R ../../standalone/Blank/gui/ standalone/${smallName}/gui
rm -R standalone/${smallName}/gui/build &> /dev/null
replace "standalone/${smallName}/gui/CMakeLists.txt"
replace "standalone/${smallName}/gui/main.cpp"
rm standalone/${smallName}/gui/core &> /dev/null
ln -sf ../cmd/core standalone/${smallName}/gui/core
ln -sf ../cmd/sandbox standalone/${smallName}/gui/sandbox

# 4) duplicate the blank paraview plugin
cp -R ../../paraview/Blank paraview/server/${Name}
mv paraview/server/${Name}/Blank.xml paraview/server/${Name}/${Name}.xml
ln -sf ../../../sandbox paraview/server/${Name}/sandbox
ln -sf ../../../../../core paraview/server/${Name}/core
replace "paraview/server/${Name}/${Name}.xml"
replace "paraview/server/${Name}/CMakeLists.txt"

echo "Project created."
