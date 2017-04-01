#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

# hard-code your custom paraview source tree path here if desired
PARAVIEW_PATH=""

# hard-code the path to your paraview installation plugins
# typical values are /usr/lib/paraview/plugins
# paraview will automatically search for plugins there and you will not have to 
# insert them manually
PARAVIEW_PLUGIN=""

CPU=""
COPY=""

UNAME=`which uname 2> /dev/null`

if [ -z "$UNAME" ]; then
  echo "Error: Please install the 'uname' package."
  exit -1
fi

while getopts "hp:P:cC" OPTION
do 
  case $OPTION in
    h)
      echo "Usage:"
      echo "  $0"
      echo "  -P [optional path to custom ParaView build]"
      echo "  -p [optional path to custom ParaView plugin directory]"
      echo "  -c (disable cpu optimization, default: on)"
      echo "  -C (force plugin copy, default: symbolic links)"
      exit 0
      ;;

    c) CPU="-DwithCpuOptimization=OFF"
      ;;

    C) COPY="COPY"
      ;;

    p)
      PARAVIEW_PLUGIN=$OPTARG
      ;;

    P)
      PARAVIEW_PATH=$OPTARG
      ;;
  esac
done

# backward compatibility
if [ -z "${PARAVIEW_PLUGIN}" ] && [ -z "${PARAVIEW_PATH}" ]; then
  if [ "$1" != "-p" ] && [ "$1" != "-c" ] && [ "$1" != "-C" ] \
  && [ "$1"  != "-P" ]; then
    PARAVIEW_PATH=$1
  fi
fi

# build server plugins
for PLUGIN in server/*; do
  if [ -d "${PLUGIN}" ]; then
    echo -e "\n\n\nBuilding plugin ${PLUGIN}..."
    cd ${PLUGIN}
    rm -R build 2> /dev/null
    mkdir build
    cd build
    if [ -z "${PARAVIEW_PATH}" ]; then
      cmake ../ ${CPU}
    else 
      cmake ../ -DParaView_DIR=${PARAVIEW_PATH} ${CPU}
    fi
    make -j 2> err.log

    # directly link the plugin in the system plugin directory if write 
    # permissions have been provided
    LIB=""
    if [ "$(uname)" == "Darwin" ]; then
      LIB=`ls lib*.dylib`
    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
      LIB=`ls lib*.so`
    elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
      # TODO: make sure it does not f**k up
      LIB=`ls lib*.dll`
    fi
    if [ ! -z "${LIB}" ]; then
      if [ ! -z "${PARAVIEW_PLUGIN}" ]; then
        rm ${PARAVIEW_PLUGIN}/${LIB} &> /dev/null
        if [ ! -z "${COPY}" ]; then
          cp `pwd`/${LIB} ${PARAVIEW_PLUGIN}/${LIB}
        else
          ln -sf `pwd`/${LIB} ${PARAVIEW_PLUGIN}/${LIB}
        fi
      fi
    fi

    ERR=`cat err.log | grep -v Processing | grep error`
    cd ../../../
    if [ ! -z "$ERR" ]; then
      echo "Build error :("
      exit
    fi
  fi
done

# build client plugins
for PLUGIN in client/*; do
  if [ -d "${PLUGIN}" ]; then
    echo "Building plugin ${PLUGIN}..."
    cd ${PLUGIN}
    rm -R build 2> /dev/null
    mkdir build
    cd build 
    if [ -z "${PARAVIEW_PATH}" ]; then
      cmake ../ ${CPU}
    else 
      cmake ../ -DParaView_DIR=${PARAVIEW_PATH} ${CPU}
    fi
    make -j 2> err.log 

    # directly link the plugin in the system plugin directory if write 
    # permissions have been provided
    LIB=""
    if [ "$(uname)" == "Darwin" ]; then
      LIB=`ls lib*.dylib`
    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
      LIB=`ls lib*.so`
    elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
      # TODO: make sure it does not f**k up
      LIB=`ls lib*.dll`
    fi
    if [ ! -z "${LIB}" ]; then
      if [ ! -z "${PARAVIEW_PLUGIN}" ]; then
        rm ${PARAVIEW_PLUGIN}/${LIB} &> /dev/null
        if [ ! -z "${COPY}" ]; then
          cp `pwd`/${LIB} ${PARAVIEW_PLUGIN}/${LIB}
        else
          ln -sf `pwd`/${LIB} ${PARAVIEW_PLUGIN}/${LIB}
        fi
      fi
    fi

    ERR=`cat err.log | grep -v Processing | grep error`
    cd ../../../
    if [ ! -z "$ERR" ]; then
      echo "Build error :("
      exit
    fi
  fi
done

echo "Success :)"
