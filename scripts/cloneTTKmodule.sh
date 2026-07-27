#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

# example: ScalarFieldSmoother MyNewModule
NameSource=$1
NameTarget=$2

smallNameTarget="$(tr '[:upper:]' '[:lower:]'  <<< "${NameTarget:0:1}")${NameTarget:1}"
bigNameTarget=$(echo "$NameTarget" | tr '[:lower:]' '[:upper:]')


smallNameSource="$(tr '[:upper:]' '[:lower:]'  <<< "${NameSource:0:1}")${NameSource:1}"
bigNameSource=$(echo "$NameSource" | tr '[:lower:]' '[:upper:]')

if [ -z "${NameTarget}" ]; then
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

echo "Cloning TTK module '${NameSource}' into '${NameTarget}'..."

function replace {
  $SED "s/${smallNameSource}/${smallNameTarget}/g" "$1" > tmp.cmake
  mv tmp.cmake "$1"
  $SED "s/${NameSource}/${NameTarget}/g" "$1" > tmp.cmake
  mv tmp.cmake "$1"
  $SED "s/${bigNameSource}/${bigNameTarget}/g" "$1" > tmp.cmake
  mv tmp.cmake "$1"
}

mkdir -p core/base 2> /dev/null
mkdir -p core/vtk 2> /dev/null
mkdir -p standalone 2> /dev/null
mkdir -p paraview/ 2> /dev/null

# 1) duplicate the source base code
echo "Creating base code functor 'core/base/${smallNameTarget}'..."
cp -R "core/base/${smallNameSource}" \
      "core/base/${smallNameTarget}"
replace "core/base/${smallNameTarget}/CMakeLists.txt"
mv "core/base/${smallNameTarget}/${NameSource}.cpp" \
   "core/base/${smallNameTarget}/${NameTarget}.cpp"
replace "core/base/${smallNameTarget}/${NameTarget}.cpp"
mv "core/base/${smallNameTarget}/${NameSource}.h" \
   "core/base/${smallNameTarget}/${NameTarget}.h"
replace "core/base/${smallNameTarget}/${NameTarget}.h"

# 2) duplicate the source wrapper
echo "Creating VTK wrapper 'core/vtk/ttk${NameTarget}'..."
cp -R "core/vtk/ttk${NameSource}" \
      "core/vtk/ttk${NameTarget}"
replace "core/vtk/ttk${NameTarget}/vtk.module"
replace "core/vtk/ttk${NameTarget}/ttk.module"
mv "core/vtk/ttk${NameTarget}/ttk${NameSource}.cpp" \
   "core/vtk/ttk${NameTarget}/ttk${NameTarget}.cpp"
replace "core/vtk/ttk${NameTarget}/ttk${NameTarget}.cpp"
mv "core/vtk/ttk${NameTarget}/ttk${NameSource}.h" \
   "core/vtk/ttk${NameTarget}/ttk${NameTarget}.h"
replace "core/vtk/ttk${NameTarget}/ttk${NameTarget}.h"

# 3) duplicate the source standalone modules
if [ -d "standalone/${NameSource}/" ]; then
  if [ ! -d "standalone/${NameSource}/cmd" ] \
    && [ ! -d "standalone/${NameSource}/gui" ]; then
    # new standalone case
    echo "Creating command line standalone program 'standalone/${NameTarget}'..."
    cp -R "standalone/${NameSource}/" \
          "standalone/${NameTarget}/"
    replace "standalone/${NameTarget}/CMakeLists.txt"
    replace "standalone/${NameTarget}/main.cpp"
  fi
  if [ -d "standalone/${NameSource}/cmd" ]; then
    echo "Creating command line standalone program 'standalone/${NameTarget}/cmd'..."
    mkdir -p "standalone/${NameTarget}"
    cp -R "standalone/${NameSource}/cmd/" \
          "standalone/${NameTarget}/cmd"
    replace "standalone/${NameTarget}/cmd/CMakeLists.txt"
    replace "standalone/${NameTarget}/cmd/main.cpp"
  fi
  if [ -d "standalone/${NameSource}/gui" ]; then
    echo "Creating GUI standalone program 'standalone/${NameTarget}/gui'..."
    mkdir -p "standalone/${NameTarget}"
    cp -R "standalone/${NameSource}/gui/" \
          "standalone/${NameTarget}/gui"
    replace "standalone/${NameTarget}/gui/CMakeLists.txt"
    replace "standalone/${NameTarget}/gui/main.cpp"
  fi
fi

# 4) duplicate the source paraview filter
echo "Creating ParaView filter 'paraview/${NameTarget}'..."
cp -R "paraview/xmls/${NameSource}.xml" \
      "paraview/xmls/${NameTarget}.xml"
replace "paraview/xmls/${NameTarget}.xml"

echo "Module ${NameTarget} created."
# Due to the use of file(GLOB...) make does not see any changes
# if cmake is not run after the new module is created
echo "/!\\ Please re-run cmake in your build folder to compile it."
