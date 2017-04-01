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
echo "Deleting project ${Name}..."

rm -R sandbox/baseCode/${smallName} 
rm -R sandbox/vtkWrappers/vtk${Name} 
rm -R standalone/${smallName} 
rm -R paraview/server/${Name} 

echo "Project deleted."
