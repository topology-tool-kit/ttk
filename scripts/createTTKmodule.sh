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
  echo "  $0 HelloWorld"
  exit -1
fi

# check for dependencies
SED=`which sed 2> /dev/null`
if [ -z "$SED" ]; then
  echo "Error: Please install sed."
  exit -1
fi

# check for paths
if [ ! -e "scripts/createTTKmodule.sh" ]; then
  echo "Error: Please run this script from the top of the source tree"
  exit -1
fi

echo "Creating TTK module ${Name}..."

function replace {
  $SED "s/blank/${smallName}/g" $1 > tmp.cmake
  mv tmp.cmake $1
  $SED "s/Blank/${Name}/g" $1 > tmp.cmake
  mv tmp.cmake $1
  $SED "s/BLANK/${bigName}/g" $1 > tmp.cmake
  mv tmp.cmake $1
}

mkdir -p core/baseCode 2> /dev/null
mkdir -p core/vtkWrappers 2> /dev/null
mkdir -p standalone 2> /dev/null
mkdir -p paraview/ 2> /dev/null

# 1) duplicate the blank base code
echo "Creating base code functor 'core/baseCode/${smallName}'"
cp -R core/baseCode/blank core/baseCode/${smallName}
mv core/baseCode/${smallName}/blank.cmake \
  core/baseCode/${smallName}/${smallName}.cmake
replace "core/baseCode/${smallName}/${smallName}.cmake"
mv core/baseCode/${smallName}/Blank.cpp \
  core/baseCode/${smallName}/${Name}.cpp
replace "core/baseCode/${smallName}/${Name}.cpp"
mv core/baseCode/${smallName}/Blank.h \
  core/baseCode/${smallName}/${Name}.h
replace "core/baseCode/${smallName}/${Name}.h"

# 2) duplicate the blank wrapper
echo "Creating VTK wrapper 'core/vtkWrappers/ttk${Name}'..."
cp -R core/vtkWrappers/ttkBlank core/vtkWrappers/ttk${Name}
mv core/vtkWrappers/ttk${Name}/ttkBlank.cmake \
  core/vtkWrappers/ttk${Name}/ttk${Name}.cmake
replace "core/vtkWrappers/ttk${Name}/ttk${Name}.cmake"
mv core/vtkWrappers/ttk${Name}/ttkBlank.cpp \
  core/vtkWrappers/ttk${Name}/ttk${Name}.cpp
replace "core/vtkWrappers/ttk${Name}/ttk${Name}.cpp"
mv core/vtkWrappers/ttk${Name}/ttkBlank.h \
  core/vtkWrappers/ttk${Name}/ttk${Name}.h
replace "core/vtkWrappers/ttk${Name}/ttk${Name}.h"

# 3) duplicate the blank standalone modules
echo "Creating command line standalone program 'standalone/${Name}/cmd'..."
mkdir standalone/${Name}
cp -R standalone/Blank/cmd/ standalone/${Name}/cmd
rm -R standalone/${Name}/cmd/build &> /dev/null
replace "standalone/${Name}/cmd/CMakeLists.txt"
replace "standalone/${Name}/cmd/main.cpp"
rm standalone/${Name}/cmd/core &> /dev/null
ln -sf ../../../core standalone/${Name}/cmd/core

echo "Creating GUI standalone program 'standalone/${Name}/gui'..."
cp -R standalone/Blank/gui/ standalone/${Name}/gui
rm -R standalone/${Name}/gui/build &> /dev/null
replace "standalone/${Name}/gui/CMakeLists.txt"
replace "standalone/${Name}/gui/main.cpp"
rm standalone/${Name}/gui/core &> /dev/null
ln -sf ..//cmd/core standalone/${Name}/gui/core

# 4) duplicate the blank paraview plugin
echo "Creating ParaView plugin 'paraview/${Name}'..."
cp -R paraview/Blank paraview/${Name}
mv paraview/${Name}/Blank.xml paraview/${Name}/${Name}.xml
ln -sf ../../core paraview/${Name}/core
replace "paraview/${Name}/${Name}.xml"
replace "paraview/${Name}/CMakeLists.txt"

# 5) update CMakeLists.txt
echo "Updating CMakeLists.txt..."
cat CMakeLists.txt | grep standalone | grep cmd > standalone-cmd.tmp
cat CMakeLists.txt | grep standalone | grep gui > standalone-gui.tmp
cat CMakeLists.txt | grep paraview > paraview.tmp

echo "cmake_minimum_required(VERSION 2.4)" > CMakeLists.txt.new
cat standalone-cmd.tmp >> CMakeLists.txt.new
echo "add_subdirectory(standalone/${Name}/cmd/)" >> CMakeLists.txt.new
cat standalone-gui.tmp >> CMakeLists.txt.new
echo "add_subdirectory(standalone/${Name}/gui/)" >> CMakeLists.txt.new
cat paraview.tmp >> CMakeLists.txt.new
echo "add_subdirectory(paraview/${Name}/)" >> CMakeLists.txt.new

mv CMakeLists.txt CMakeLists.txt.old
mv CMakeLists.txt.new CMakeLists.txt
rm standalone-cmd.tmp 2> /dev/null
rm standalone-gui.tmp 2> /dev/null
rm paraview.tmp 2> /dev/null

echo "Module created."
