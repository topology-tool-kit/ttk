#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

# example: ScalarFieldSmoother
Name=$1

if [ -z "${Name}" ]; then
  echo "Usage:"
  echo "  $0 <Name, first letter in uppercase, no space>"
  echo "Example:"
  echo "  $0 MyNewModule"
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

echo "Creating TTK module ${Name}..."

scripts/cloneTTKmodule.sh HelloWorld "${Name}"
