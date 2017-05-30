#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

function print_usage(){

  echo "Usage:"
  echo "  $0"
  echo "  -n <Name of the module, first letter uppercase, no space>"
  echo "  -a [Optional code anonymization]"
  echo "  -c [With command line program]"
  echo "  -d [With documentation]"
  echo "  -g [With graphical user interface program]"
  echo "  -p [With ParaView plugin]"

  exit 0
}

function anonymize(){

  if [ -d $1 ]; then
    # quit if directory
    return 0;
  fi

  echo "Anonymizing file '$1'..."
  cat $1 > tmp2
  cat tmp2 | grep --binary-files=text -v "\author" > tmp
  cat tmp | grep --binary-files=text -v "author:" > tmp2
  cat tmp2 | grep --binary-files=text -v "Copyright" > tmp
  cat tmp | grep --binary-files=text -v "copyright" > tmp2
  mv tmp2 $1
  rm tmp
}

function process_dependency(){

  local CMAKE=$1
  local EXPORT="1"

  if [ -z `echo $1 | grep CMakeLists.txt` ]; then
    # not the root, change the filename
    DIR="$1"
    if [ "$DIR" == "core" ]; then
      DIR=$CORE_PATH
    fi
    CMAKE="$DIR/$2/$3/$3.cmake"

    # export the actual entry
    if [ ! -e release/$MODULE/ttk/$2/$3 ]; then
      echo "Processing $1 $2 dependency '$3'..."
      if [ "$1" == "core" ] && [ "$2" == "baseCode" ]; then
        cp -R $CORE_PATH/baseCode/$3 \
          release/$MODULE/ttk/$2/$3 > /dev/null
      fi
      if [ "$1" == "core" ] && [ "$2" == "vtkWrappers" ]; then
        cp -R $CORE_PATH/vtkWrappers/$3 \
          release/$MODULE/ttk/$2/$3 > /dev/null
      fi
    else
      EXPORT=""
    fi
  else
    echo "Processing dependencies for module's root: $1"
  fi

  if [ ! -z "$EXPORT" ]; then

    # baseCode dependencies
    for i in `cat $CMAKE | grep baseCode_package`; do
      entry=`echo $i | cut -d '(' -f2 | cut -d ')' -f1`
      if [ -e $CORE_PATH/baseCode/$entry ]; then
        process_dependency \
          core \
          baseCode \
          $entry
      else
        echo "Error! Could not find dependency '$entry'!"
        exit -1
      fi
    done

    # vtkWrapper dependencies
    for i in `cat $CMAKE | grep vtkWrapper_package`; do
      entry=`echo $i | cut -d '(' -f2 | cut -d ')' -f1`
      if [ -e $CORE_PATH/vtkWrappers/$entry ]; then
        process_dependency \
          core \
          vtkWrappers \
          $entry
      else
        echo "Error! Could not find dependency '$entry'!"
        exit -1
      fi
    done
  fi

  if [ `echo $1 | grep CMakeLists.txt` ]; then
    # finishing the dependencies for the root
    cp -R $CORE_PATH/ttk.cmake \
      release/$MODULE/ttk/ttk.cmake > /dev/null
    cp -R $CORE_PATH/ttk.doxygen \
      release/$MODULE/ttk/ttk.doxygen > /dev/null
    cp -R $CORE_PATH/baseCode/baseCode.cmake \
      release/$MODULE/ttk/baseCode/baseCode.cmake > /dev/null

    cp -R $CORE_PATH/baseCode/common \
      release/$MODULE/ttk/baseCode/common > /dev/null
    cp -R $CORE_PATH/baseCode/triangulation \
      release/$MODULE/ttk/baseCode/triangulation > /dev/null
  fi
}

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


ANONYMOUS=""
COMMANDLINE=""
GUI=""
PARAVIEW=""
MODULE=""
ROOT=""
DOCUMENTATION=""
IS_SANDBOX=""
CORE_PATH="core/"
LICENSE="LICENSE"

if [ ! -e "$CORE_PATH" ]; then
  CORE_PATH="../../core/"
  IS_SANDBOX="1"
fi

while getopts "acdghn:p" OPTION
do
  case $OPTION in
    a)
      ANONYMOUS="1"
      ;;
    c)
      COMMANDLINE="1"
      ;;
    d)
      DOCUMENTATION="1"
      ;;
    g)
      GUI="1"
      ;;
    h)
      print_usage
      ;;
    n)
      MODULE=$OPTARG
      ;;
    p)
      PARAVIEW="1"
      ;;
  esac 
done

if [ -z "$MODULE" ]; then
  print_usage
fi

if [ -z "$COMMANDLINE" ] && [ -z "$PARAVIEW" ] && [ -z "$GUI" ]; then
  print_usage
fi

smallName=${MODULE,}
Name=${MODULE}
bigName=${Name^^}

if [ ! -z "$PARAVIEW" ]; then
  ROOT=paraview/$MODULE/
fi

if [ ! -z "$COMMANDLINE" ]; then
  ROOT=standalone/$MODULE/cmd/
fi

if [ ! -z "$GUI" ]; then
  ROOT=standalone/$MODULE/gui/
fi

if [ ! -e $ROOT/CMakeLists.txt ]; then
  echo "Error: Could not find module's root '${ROOT}/CMakeLists.txt'!"
  exit -1
fi

if [ -e release/$MODULE ]; then
  echo "Removing previous instance of module release..."
  rm -R release/$MODULE &> /dev/null
fi 

# 1. create the main folder
mkdir -p release/$MODULE/data
mkdir -p release/$MODULE/ttk/baseCode
mkdir -p release/$MODULE/ttk/vtkWrappers

# 2. Process the code dependencies
process_dependency $ROOT/CMakeLists.txt

# 3. Create the paraview plugin
if [ ! -z "$PARAVIEW" ]; then
  echo "Creating ParaView plugin..."
  mkdir release/$MODULE/paraview
  cp paraview/$MODULE/CMakeLists.txt release/$MODULE/paraview
  cp paraview/$MODULE/*xml release/$MODULE/paraview
  $SED "s/core/ttk/g" release/$MODULE/paraview/CMakeLists.txt > tmp.txt
  mv tmp.txt release/$MODULE/paraview/CMakeLists.txt
  echo "  Assuming binaries (if any) have been compiled in a 'build' subdir..."
  rm -R release/$MODULE/paraview/build &> /dev/null
fi

# 4. Create the command line program
if [ ! -z "$COMMANDLINE" ]; then
  echo "Creating command line program..."
  cp -R $CORE_PATH/vtkWrappers/ttkWrapper \
    release/${MODULE}/ttk/vtkWrappers/ttkWrapper > /dev/null
  cp -R $CORE_PATH/vtkWrappers/ttkTriangulation \
    release/${MODULE}/ttk/vtkWrappers/ttkTriangulation > /dev/null
  cp -R $CORE_PATH/vtkWrappers/ttkProgramBase \
    release/$MODULE/ttk/vtkWrappers/ttkProgramBase > /dev/null
  cp -R standalone/$Name/cmd release/$MODULE/cmd > /dev/null
  $SED "s/core/ttk/g" release/$MODULE/cmd/CMakeLists.txt > tmp.txt
  mv tmp.txt release/$MODULE/cmd/CMakeLists.txt
  rm release/$MODULE/cmd/core &> /dev/null
  rm release/$MODULE/cmd/data &> /dev/null
  echo "  Assuming binaries (if any) have been compiled in a 'build' subdir..."
  rm -R release/$MODULE/cmd/build &> /dev/null
fi

# 5. Create the gui program
if [ ! -z "$GUI" ]; then
  echo "Creating graphical user interface program..."
  cp -R $CORE_PATH/vtkWrappers/ttkWrapper \
    release/${MODULE}/ttk/vtkWrappers/ttkWrapper > /dev/null
  cp -R $CORE_PATH/vtkWrappers/ttkTriangulation \
    release/${MODULE}/ttk/vtkWrappers/ttkTriangulation > /dev/null
  cp -R $CORE_PATH/vtkWrappers/ttkUserInterfaceBase \
    release/$MODULE/ttk/vtkWrappers/ttkUserInterfaceBase > /dev/null
  cp -R $CORE_PATH/vtkWrappers/ttkTextureMapFromField \
    release/$MODULE/ttk/vtkWrappers/ttkTextureMapFromField > /dev/null
  cp -R $CORE_PATH/vtkWrappers/ttkWRLExporter \
    release/$MODULE/ttk/vtkWrappers/ttkWRLExporter > /dev/null
  cp -R standalone/$Name/gui release/$MODULE/gui > /dev/null
  $SED "s/core/ttk/g" release/$MODULE/gui/CMakeLists.txt > tmp.txt
  mv tmp.txt release/$MODULE/gui/CMakeLists.txt
  rm release/$MODULE/gui/core
  rm release/$MODULE/gui/data &> /dev/null
  rm release/$MODULE/gui/textures &> /dev/null
  mkdir -p release/$MODULE/gui/textures
  cp -R standalone/$Name/gui/textures/png \
    release/$MODULE/gui/textures
  echo "  Assuming binaries (if any) have been compiled in a 'build' subdir..."
  rm -R release/$MODULE/gui/build &> /dev/null
fi

# 6. Create license and readme
cp $LICENSE release/$MODULE/LICENSE

# 7. Anonymize code
if [ ! -z "$ANONYMOUS" ]; then
  if [ ! -z "$COMMANDLINE" ]; then
    for i in release/$MODULE/cmd/*; do
      anonymize $i
    done
  fi
  if [ ! -z "$GUI" ]; then
    for i in release/$MODULE/gui/*; do
      anonymize $i
    done
  fi
  if [ ! -z "$PARAVIEW" ]; then
    for i in release/$MODULE/paraview/*; do
      anonymize $i
    done
  fi
  for i in release/$MODULE/ttk/*; do
    anonymize $i
  done
  for i in release/$MODULE/ttk/baseCode/*; do
    anonymize $i
  done
  for i in release/$MODULE/ttk/vtkWrappers/*; do
    anonymize $i
  done
  for i in release/$MODULE/ttk/baseCode/*/*; do
    anonymize $i
  done
  for i in release/$MODULE/ttk/vtkWrappers/*/*; do
    anonymize $i
  done
fi

# 8. Create the documentation
if [ ! -z "$DOCUMENTATION" ]; then
  echo "Generating documentation..."
  cd release/$MODULE/ttk
  doxygen ttk.doxygen &> /dev/null
  cd ../../../
  mv release/doc release/$MODULE/doc/
fi 

# Message: provide data-sets, update README file and zip and ready to ship!
echo ""
echo "Release directory 'release/$MODULE' completed!"
echo ""
echo "Now:"
echo "=="
echo "  1) Copy some carefully-chosen test data-sets to the following folder:"
echo "    release/$MODULE/data"
echo "  2) Create and edit your README file:"
echo "    release/$MODULE/README"
echo "  3) Edit the default LICENSE file:"
echo "    release/$MODULE/LICENSE"
echo "  4) Once you are finished, zip the release directory:"
echo "    release/$MODULE"
echo ""
echo "And that's it! Your release is all set to ship!"
