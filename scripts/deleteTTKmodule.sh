#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

# example: ScalarFieldSmoother
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

# check for paths
if [ ! -e "scripts/createTTKmodule.sh" ]; then
  echo "Error: Please run this script from the top of the source tree"
  exit -1
fi

# check for dependencies
echo "Deleting TTK module ${Name}..."

rm -R core/base/${smallName}  2> /dev/null
rm -R core/vtk/ttk${Name} 2> /dev/null
rm -R standalone/${Name} 2> /dev/null
rm -R paraview/${Name} 2> /dev/null

echo "Module deleted."
