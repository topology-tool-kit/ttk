#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

function print_usage(){

  echo "Usage:"
  echo "  $0"
  echo "  -n <Name of the project, first letter uppercase, no space>"
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
    if [ ! -e release/$PROJECT/ttk/$2/$3 ]; then
      echo "Processing $1 $2 dependency '$3'..."
      if [ "$1" == "core" ] && [ "$2" == "baseCode" ]; then
        cp -R $CORE_PATH/baseCode/$3 \
          release/$PROJECT/ttk/$2/$3 > /dev/null
      fi
      if [ "$1" == "core" ] && [ "$2" == "vtkWrappers" ]; then
        cp -R $CORE_PATH/vtkWrappers/$3 \
          release/$PROJECT/ttk/$2/$3 > /dev/null
      fi
      if [ "$1" == "sandbox" ] && [ "$2" == "baseCode" ]; then
        cp -R sandbox/baseCode/$3 \
          release/$PROJECT/ttk/$2/$3 > /dev/null
      fi
      if [ "$1" == "sandbox" ] && [ "$2" == "vtkWrappers" ]; then
        cp -R sandbox/vtkWrappers/$3 \
          release/$PROJECT/ttk/$2/$3 > /dev/null
      fi
      
    else
      EXPORT=""
    fi
  else
    echo "Processing dependencies for project's root: $1"
  fi

  if [ ! -z "$EXPORT" ]; then

    # baseCode dependencies
    for i in `cat $CMAKE | grep baseCode_package`; do
      entry=`echo $i | cut -d '(' -f2 | cut -d ')' -f1`
      if [ -e sandbox/baseCode/$entry/ ]; then
        process_dependency \
          sandbox \
          baseCode \
          $entry
      else
        if [ -e $CORE_PATH/baseCode/$entry ]; then
          process_dependency \
            core \
            baseCode \
            $entry
        else
          echo "Error! Could not find dependency '$entry'!"
          exit -1
        fi
      fi
    done

    # vtkWrapper dependencies
    for i in `cat $CMAKE | grep vtkWrapper_package`; do
      entry=`echo $i | cut -d '(' -f2 | cut -d ')' -f1`
      if [ -e sandbox/vtkWrappers/$entry/ ]; then
        process_dependency \
          sandbox \
          vtkWrappers \
          $entry
      else
        if [ -e $CORE_PATH/vtkWrappers/$entry ]; then
          process_dependency \
            core \
            vtkWrappers \
            $entry
        else
          echo "Error! Could not find dependency '$entry'!"
          exit -1
        fi
      fi
    done
  fi

  if [ `echo $1 | grep CMakeLists.txt` ]; then
    # finishing the dependencies for the root
    cp -R $CORE_PATH/ttk.cmake \
      release/$PROJECT/ttk/ttk.cmake > /dev/null
    cp -R $CORE_PATH/ttk.doxygen \
      release/$PROJECT/ttk/ttk.doxygen > /dev/null
    cp -R $CORE_PATH/baseCode/baseCode.cmake \
      release/$PROJECT/ttk/baseCode/baseCode.cmake > /dev/null

    cp -R $CORE_PATH/baseCode/common \
      release/$PROJECT/ttk/baseCode/common > /dev/null
    cp -R $CORE_PATH/baseCode/triangulation \
      release/$PROJECT/ttk/baseCode/triangulation > /dev/null
  fi
}

# check for dependencies
SED=`which sed 2> /dev/null`
if [ -z "$SED" ]; then
  echo "Error: Please install sed."
  exit -1
fi

ANONYMOUS=""
COMMANDLINE=""
GUI=""
PARAVIEW=""
PROJECT=""
ROOT=""
DOCUMENTATION=""
IS_SANDBOX=""
CORE_PATH="core/"
LICENSE="../../LICENSE"

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
      PROJECT=$OPTARG
      ;;
    p)
      PARAVIEW="1"
      ;;
  esac 
done

if [ -z "$PROJECT" ]; then
  print_usage
fi

if [ -z "$COMMANDLINE" ] && [ -z "$PARAVIEW" ] && [ -z "$GUI" ]; then
  print_usage
fi

smallName=${PROJECT,}
bigName=${Name^^}



# assuming server plugins
if [ ! -z "$PARAVIEW" ]; then
  ROOT=paraview/server/$PROJECT/
fi

if [ ! -z "$COMMANDLINE" ]; then
  ROOT=standalone/$smallName/cmd/
fi

if [ ! -z "$GUI" ]; then
  ROOT=standalone/$smallName/gui/
fi

if [ ! -e $ROOT/CMakeLists.txt ]; then
  echo "Error: Could not find project's root '${ROOT}/CMakeLists.txt'!"
  exit -1
fi

if [ -e release/$PROJECT ]; then
  echo "Removing previous instance of project release..."
  rm -R release/$PROJECT &> /dev/null
fi 

# 1. create the main folder
mkdir -p release/$PROJECT/data
mkdir -p release/$PROJECT/ttk/baseCode
mkdir -p release/$PROJECT/ttk/vtkWrappers

LICENSE="../../LICENSE"

# 2. Process the code dependencies
process_dependency $ROOT/CMakeLists.txt

# 3. Create the paraview plugin
if [ ! -z "$PARAVIEW" ]; then
  # assuming server plugins
  echo "Creating ParaView plugin..."
  mkdir release/$PROJECT/paraview
  cp paraview/server/$PROJECT/CMakeLists.txt release/$PROJECT/paraview
  cp paraview/server/$PROJECT/*xml release/$PROJECT/paraview
  $SED "s/core/ttk/g" release/$PROJECT/paraview/CMakeLists.txt > tmp.txt
  mv tmp.txt release/$PROJECT/paraview/CMakeLists.txt
  echo "  Assuming binaries (if any) have been compiled in a 'build' subdir..."
  rm -R release/$PROJECT/paraview/build &> /dev/null
fi

# 4. Create the command line program
if [ ! -z "$COMMANDLINE" ]; then
  # assuming server plugins
  echo "Creating command line program..."
  cp -R $CORE_PATH/vtkWrappers/ttkWrapper \
    release/${PROJECT}/ttk/vtkWrappers/ttkWrapper > /dev/null
  cp -R $CORE_PATH/vtkWrappers/vtkTriangulation \
    release/${PROJECT}/ttk/vtkWrappers/vtkTriangulation > /dev/null
  cp -R $CORE_PATH/vtkWrappers/vtkProgramBase \
    release/$PROJECT/ttk/vtkWrappers/vtkProgramBase > /dev/null
  cp -R standalone/$smallName/cmd release/$PROJECT/cmd > /dev/null
  $SED "s/core/ttk/g" release/$PROJECT/cmd/CMakeLists.txt > tmp.txt
  mv tmp.txt release/$PROJECT/cmd/CMakeLists.txt
  rm release/$PROJECT/cmd/core &> /dev/null
  rm release/$PROJECT/cmd/data &> /dev/null
  rm release/$PROJECT/cmd/sandbox &> /dev/null
  echo "  Assuming binaries (if any) have been compiled in a 'build' subdir..."
  rm -R release/$PROJECT/cmd/build &> /dev/null
fi

# 5. Create the gui program
if [ ! -z "$GUI" ]; then
  # assuming server plugins
  echo "Creating graphical user interface program..."
  cp -R $CORE_PATH/vtkWrappers/ttkWrapper \
    release/${PROJECT}/ttk/vtkWrappers/ttkWrapper > /dev/null
  cp -R $CORE_PATH/vtkWrappers/vtkTriangulation \
    release/${PROJECT}/ttk/vtkWrappers/vtkTriangulation > /dev/null
  cp -R $CORE_PATH/vtkWrappers/vtkUserInterfaceBase \
    release/$PROJECT/ttk/vtkWrappers/vtkUserInterfaceBase > /dev/null
  cp -R $CORE_PATH/vtkWrappers/vtkTextureMapFromField \
    release/$PROJECT/ttk/vtkWrappers/vtkTextureMapFromField > /dev/null
  cp -R $CORE_PATH/vtkWrappers/vtkWRLExporter \
    release/$PROJECT/ttk/vtkWrappers/vtkWRLExporter > /dev/null
  cp -R standalone/$smallName/gui release/$PROJECT/gui > /dev/null
  $SED "s/core/ttk/g" release/$PROJECT/gui/CMakeLists.txt > tmp.txt
  mv tmp.txt release/$PROJECT/gui/CMakeLists.txt
  rm release/$PROJECT/gui/core
  rm release/$PROJECT/gui/data &> /dev/null
  rm release/$PROJECT/gui/sandbox &> /dev/null
  rm release/$PROJECT/gui/textures &> /dev/null
  mkdir -p release/$PROJECT/gui/textures
  cp -R standalone/$smallName/gui/textures/png \
    release/$PROJECT/gui/textures
  echo "  Assuming binaries (if any) have been compiled in a 'build' subdir..."
  rm -R release/$PROJECT/gui/build &> /dev/null
fi

# 6. Create license and readme
cp $LICENSE release/$PROJECT/LICENSE

# 7. Anonymize code
if [ ! -z "$ANONYMOUS" ]; then
  if [ ! -z "$COMMANDLINE" ]; then
    for i in release/$PROJECT/cmd/*; do
      anonymize $i
    done
  fi
  if [ ! -z "$GUI" ]; then
    for i in release/$PROJECT/gui/*; do
      anonymize $i
    done
  fi
  if [ ! -z "$PARAVIEW" ]; then
    for i in release/$PROJECT/paraview/*; do
      anonymize $i
    done
  fi
  for i in release/$PROJECT/ttk/*; do
    anonymize $i
  done
  for i in release/$PROJECT/ttk/baseCode/*; do
    anonymize $i
  done
  for i in release/$PROJECT/ttk/vtkWrappers/*; do
    anonymize $i
  done
  for i in release/$PROJECT/ttk/baseCode/*/*; do
    anonymize $i
  done
  for i in release/$PROJECT/ttk/vtkWrappers/*/*; do
    anonymize $i
  done
fi

# 8. Create the documentation
if [ ! -z "$DOCUMENTATION" ]; then
  echo "Generating documentation..."
  cd release/$PROJECT/ttk
  doxygen ttk.doxygen &> /dev/null
  cd ../../../
  mv release/doc release/$PROJECT/doc/
fi 

# Message: provide data-sets, update README file and zip and ready to ship!
echo ""
echo "Release directory 'release/$PROJECT' completed!"
echo ""
echo "Now:"
echo "=="
echo "  1) Copy some carefully-chosen test data-sets to the following folder:"
echo "    release/$PROJECT/data"
echo "  2) Create and edit your README file:"
echo "    release/$PROJECT/README"
echo "  3) Edit the default LICENSE file:"
echo "    release/$PROJECT/LICENSE"
echo "  4) Once you are finished, zip the release directory:"
echo "    release/$PROJECT"
echo ""
echo "And that's it! Your release is all set to ship!"
