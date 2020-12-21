#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

# example: ScalarFieldSmoother MyNewModule
NameSource=$1
NameDestination=$2

smallNameDestination="$(tr '[:upper:]' '[:lower:]'  <<< "${NameDestination:0:1}")${NameDestination:1}"
bigNameDestination=$(echo "$NameDestination" | tr '[:lower:]' '[:upper:]')


smallNameSource="$(tr '[:upper:]' '[:lower:]'  <<< "${NameSource:0:1}")${NameSource:1}"
bigNameSource=$(echo "$NameSource" | tr '[:lower:]' '[:upper:]')

if [ -z "${NameDestination}" ]; then
  echo "Usage:"
  echo "  $0 <SourceTTKmodule> <DestinationTTKmodule>"
  echo "Example:"
  echo "  $0 ScalarFieldSmoother MyNewModule"
  exit 1
fi

# check for dependencies
SED=$(command -v sed 2> /dev/null)
if [ -z "$SED" ]; then
  echo "Error: Please install sed."
  exit 2
fi

# check for paths
if [ ! -e "scripts/createTTKmodule.sh" ]; then
  echo "Error: Please run this script from the top of the source tree"
  exit 3
fi

echo "Cloning TTK module '${NameSource}' into '${NameDestination}'..."

function replace {
  $SED "s/${smallNameSource}/${smallNameDestination}/g" "$1" > tmp.cmake
  mv tmp.cmake "$1"
  $SED "s/${NameSource}/${NameDestination}/g" "$1" > tmp.cmake
  mv tmp.cmake "$1"
  $SED "s/${bigNameSource}/${bigNameDestination}/g" "$1" > tmp.cmake
  mv tmp.cmake "$1"
}

mkdir -p core/base 2> /dev/null
mkdir -p core/vtk 2> /dev/null
mkdir -p standalone 2> /dev/null
mkdir -p paraview/ 2> /dev/null

# 1) duplicate the source base code
echo "Creating base code functor 'core/base/${smallNameDestination}'..."
cp -R "core/base/${smallNameSource}" \
      "core/base/${smallNameDestination}"
replace "core/base/${smallNameDestination}/CMakeLists.txt"
mv "core/base/${smallNameDestination}/${NameSource}.cpp" \
   "core/base/${smallNameDestination}/${NameDestination}.cpp"
replace "core/base/${smallNameDestination}/${NameDestination}.cpp"
mv "core/base/${smallNameDestination}/${NameSource}.h" \
   "core/base/${smallNameDestination}/${NameDestination}.h"
replace "core/base/${smallNameDestination}/${NameDestination}.h"

# 2) duplicate the source wrapper
echo "Creating VTK wrapper 'core/vtk/ttk${NameDestination}'..."
cp -R "core/vtk/ttk${NameSource}" \
      "core/vtk/ttk${NameDestination}"
replace "core/vtk/ttk${NameDestination}/vtk.module"
replace "core/vtk/ttk${NameDestination}/ttk.module"
mv "core/vtk/ttk${NameDestination}/ttk${NameSource}.cpp" \
   "core/vtk/ttk${NameDestination}/ttk${NameDestination}.cpp"
replace "core/vtk/ttk${NameDestination}/ttk${NameDestination}.cpp"
mv "core/vtk/ttk${NameDestination}/ttk${NameSource}.h" \
   "core/vtk/ttk${NameDestination}/ttk${NameDestination}.h"
replace "core/vtk/ttk${NameDestination}/ttk${NameDestination}.h"

# 3) duplicate the source standalone modules
if [ -d "standalone/${NameSource}/" ]; then
  if [ ! -d "standalone/${NameSource}/cmd" ] \
    && [ ! -d "standalone/${NameSource}/gui" ]; then
    # new standalone case
    echo "Creating command line standalone program 'standalone/${NameDestination}'..."
    cp -R "standalone/${NameSource}/" \
          "standalone/${NameDestination}/"
    replace "standalone/${NameDestination}/CMakeLists.txt"
    replace "standalone/${NameDestination}/main.cpp"
  fi
  if [ -d "standalone/${NameSource}/cmd" ]; then
    echo "Creating command line standalone program 'standalone/${NameDestination}/cmd'..."
    mkdir -p "standalone/${NameDestination}"
    cp -R "standalone/${NameSource}/cmd/" \
          "standalone/${NameDestination}/cmd"
    replace "standalone/${NameDestination}/cmd/CMakeLists.txt"
    replace "standalone/${NameDestination}/cmd/main.cpp"
  fi
  if [ -d "standalone/${NameSource}/gui" ]; then
    echo "Creating GUI standalone program 'standalone/${NameDestination}/gui'..."
    mkdir -p "standalone/${NameDestination}"
    cp -R "standalone/${NameSource}/gui/" \
          "standalone/${NameDestination}/gui"
    replace "standalone/${NameDestination}/gui/CMakeLists.txt"
    replace "standalone/${NameDestination}/gui/main.cpp"
  fi
fi

# 4) duplicate the source paraview filter
echo "Creating ParaView filter 'paraview/${NameDestination}'..."
cp -R "paraview/${NameSource}/" \
      "paraview/${NameDestination}"
replace "paraview/${NameDestination}/TTKFilter.cmake"
mv "paraview/${NameDestination}/${NameSource}.xml" \
   "paraview/${NameDestination}/${NameDestination}.xml"
replace "paraview/${NameDestination}/${NameDestination}.xml"

echo "Module ${NameDestination} created."
# Due to the use of file(GLOB...) make does not see any changes
# if cmake is not run after the new module is created
echo "/!\\ Please re-run cmake in your build folder to compile it."
