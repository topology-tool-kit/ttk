#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

# example: ScalarFieldSmoother HelloWorld
NameSource=$1
NameDestination=$2

smallNameDestination="`tr '[:upper:]' '[:lower:]'  <<< ${NameDestination:0:1}`${NameDestination:1}"
bigNameDestination=`echo $NameDestination | tr '[:lower:]' '[:upper:]'`


smallNameSource="`tr '[:upper:]' '[:lower:]'  <<< ${NameSource:0:1}`${NameSource:1}"
bigNameSource=`echo $NameSource | tr '[:lower:]' '[:upper:]'`

if [ -z "${NameDestination}" ]; then
  echo "Usage:"
  echo "  $0 <SourceTTKmodule> <DestinationTTKmodule>"
  echo "Example:"
  echo "  $0 ScalarFieldSmoother HelloWorld"
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

echo "Cloning TTK module '${NameSource}' into '${NameDestination}'..."

function replace {
  $SED "s/${smallNameSource}/${smallNameDestination}/g" $1 > tmp.cmake
  mv tmp.cmake $1
  $SED "s/${NameSource}/${NameDestination}/g" $1 > tmp.cmake
  mv tmp.cmake $1
  $SED "s/${bigNameSource}/${bigNameDestination}/g" $1 > tmp.cmake
  mv tmp.cmake $1
}

mkdir -p core/base 2> /dev/null
mkdir -p core/vtk 2> /dev/null
mkdir -p standalone 2> /dev/null
mkdir -p paraview/ 2> /dev/null

# 1) duplicate the blank base code
echo "Creating base code functor 'core/base/${smallNameDestination}'"
cp -R core/base/${smallNameSource} core/base/${smallNameDestination}
replace "core/base/${smallNameDestination}/CMakeLists.txt"
mv core/base/${smallNameDestination}/${NameSource}.cpp \
  core/base/${smallNameDestination}/${NameDestination}.cpp
replace "core/base/${smallNameDestination}/${NameDestination}.cpp"
mv core/base/${smallNameDestination}/${NameSource}.h \
  core/base/${smallNameDestination}/${NameDestination}.h
replace "core/base/${smallNameDestination}/${NameDestination}.h"

# 2) duplicate the blank wrapper
echo "Creating VTK wrapper 'core/vtk/ttk${NameDestination}'..."
cp -R core/vtk/ttk${NameSource} core/vtk/ttk${NameDestination}
replace "core/vtk/ttk${NameDestination}/CMakeLists.txt"
mv core/vtk/ttk${NameDestination}/ttk${NameSource}.cpp \
  core/vtk/ttk${NameDestination}/ttk${NameDestination}.cpp
replace "core/vtk/ttk${NameDestination}/ttk${NameDestination}.cpp"
mv core/vtk/ttk${NameDestination}/ttk${NameSource}.h \
  core/vtk/ttk${NameDestination}/ttk${NameDestination}.h
replace "core/vtk/ttk${NameDestination}/ttk${NameDestination}.h"

# 3) duplicate the blank standalone modules
echo "Creating command line standalone program 'standalone/${NameDestination}/cmd'..."
mkdir standalone/${NameDestination}
cp -R standalone/${NameSource}/cmd/ standalone/${NameDestination}/cmd
rm -R standalone/${NameDestination}/cmd/build &> /dev/null
replace "standalone/${NameDestination}/cmd/CMakeLists.txt"
replace "standalone/${NameDestination}/cmd/main.cpp"

echo "Creating GUI standalone program 'standalone/${NameDestination}/gui'..."
cp -R standalone/${NameSource}/gui/ standalone/${NameDestination}/gui
rm -R standalone/${NameDestination}/gui/build &> /dev/null
replace "standalone/${NameDestination}/gui/CMakeLists.txt"
replace "standalone/${NameDestination}/gui/main.cpp"

# 4) duplicate the blank paraview plugin
echo "Creating ParaView plugin 'paraview/${NameDestination}'..."
cp -R paraview/${NameSource}/ paraview/${NameDestination}
mv paraview/${NameDestination}/${NameSource}.xml paraview/${NameDestination}/${NameDestination}.xml
replace "paraview/${NameDestination}/${NameDestination}.xml"
replace "paraview/${NameDestination}/CMakeLists.txt"


echo "Module created."
