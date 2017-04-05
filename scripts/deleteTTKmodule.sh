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

# check for paths
if [ ! -e "scripts/createTTKmodule.sh" ]; then
  echo "Error: Please run this script from the top of the source tree"
  exit -1
fi

# check for dependencies
echo "Deleting TTK module ${Name}..."

rm -R core/baseCode/${smallName}  2> /dev/null
rm -R core/vtkWrappers/vtk${Name} 2> /dev/null
rm -R standalone/${Name} 2> /dev/null
rm -R paraview/${Name} 2> /dev/null

cat CMakeLists.txt | grep -v ${Name} > CMakeLists.txt.new
mv CMakeLists.txt CMakeLists.txt.old
mv CMakeLists.txt.new CMakeLists.txt

echo "Module deleted."
