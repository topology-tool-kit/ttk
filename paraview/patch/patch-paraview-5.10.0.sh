#!/bin/bash
# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>
PATCH_BIN=$(command -v patch 2> /dev/null)

if [ -z "$PATCH_BIN" ]; then
  echo "Error: Please install the 'patch' package."
  exit 1
fi

if [ -z "$1" ] ||[ ! -d "$1" ]; then
  echo "Error: Missing ParaView source tree path."
  echo "Usage:"
  echo "   ./patch.sh <ParaView source tree path>"
  exit 2
fi

PATCH_DIR=$(pwd)

if [ ! -f ${PATCH_DIR}/paraview-5.8.0-CMakeLists.txt.patch ]; then
  echo "You should run this script directly from the paraview/patch folder"
  exit 3
fi

cd "$1" || exit 3
$PATCH_BIN Clients/ParaView/CMakeLists.txt \
  < "${PATCH_DIR}/paraview-5.8.0-CMakeLists.txt.patch"
$PATCH_BIN Qt/Components/Resources/UI/pqAboutDialog.ui \
  < "${PATCH_DIR}/pqAboutDialog.ui.patch"
$PATCH_BIN VTK/IO/Export/vtkVRMLExporter.cxx \
  < "${PATCH_DIR}/paraview-5.8.0-vtkVRMLExporter.cxx.patch"
$PATCH_BIN VTK/Filters/Extraction/vtkExtractSelectedIds.cxx \
  < "${PATCH_DIR}/paraview-5.8.0-vtkExtractSelectedIds.cxx.patch"

cp "${PATCH_DIR}/splash.png"  Clients/ParaView/PVSplashScreen.png
cp "${PATCH_DIR}/icon128.png" Clients/ParaView/pvIcon.png
cp "${PATCH_DIR}/icon512.png" Clients/ParaView/pvIcon-512x512.png
cp "${PATCH_DIR}/icon96.png"  Clients/ParaView/pvIcon-96x96.png
cp "${PATCH_DIR}/icon64.png"  Clients/ParaView/pvIcon-64x64.png
cp "${PATCH_DIR}/icon32.png"  Clients/ParaView/pvIcon-32x32.png
cp "${PATCH_DIR}/icon22.png"  Clients/ParaView/pvIcon-22x22.png
cp "${PATCH_DIR}/icon16.png"  Clients/ParaView/pvIcon-16x16.png
cp "${PATCH_DIR}/icon.ico"    Clients/ParaView/WinIcon.ico
cp "${PATCH_DIR}/icon.ico"    Clients/ParaView/pvIcon.ico

cp "${PATCH_DIR}/splash.png" \
  Qt/Components/Resources/Icons/PVSplashScreen.png
cp "${PATCH_DIR}/icon.ico" \
  Qt/Components/Resources/Icons/paraqlogo.ico
cp "${PATCH_DIR}/icon16.png" \
  Qt/Components/Resources/Icons/pqAppIcon16.png
cp "${PATCH_DIR}/icon22.png" \
  Qt/Components/Resources/Icons/pqAppIcon22.png
cp "${PATCH_DIR}/icon32.png" \
  Qt/Components/Resources/Icons/pqAppIcon32.png
cp "${PATCH_DIR}/icon64.png" \
  Qt/Components/Resources/Icons/pqAppIcon64.png
cp "${PATCH_DIR}/icon512.png" \
  Qt/Components/Resources/Icons/pvIcon512.png
cp "${PATCH_DIR}/icon96.png" \
  Qt/Components/Resources/Icons/pvIcon96.png
cp "${PATCH_DIR}/icon64.png" \
  Qt/Components/Resources/Icons/pvIcon64.png
cp "${PATCH_DIR}/icon32.png" \
  Qt/Components/Resources/Icons/pvIcon32.png
cp "${PATCH_DIR}/icon22.png" \
  Qt/Components/Resources/Icons/pvIcon22.png
cp "${PATCH_DIR}/icon16.png" \
  Qt/Components/Resources/Icons/pvIcon16.png

# # processing example data-sets
mkdir -p TTK/Data/
cp ${PATCH_DIR}/data/* TTK/Data/
cp TTK/Data/*pvsm Qt/ApplicationComponents/Resources/ExampleVisualizations/
cp TTK/Data/*png Qt/ApplicationComponents/Resources/Thumbnails/
$PATCH_BIN CMakeLists.txt \
  < "${PATCH_DIR}/paraview-examples-5.10-CMakeLists.txt.patch"
$PATCH_BIN Qt/ApplicationComponents/Resources/pqApplicationComponents.qrc \
  < "${PATCH_DIR}/paraview-5.5.0-pqApplicationComponents.qrc.patch"
$PATCH_BIN Qt/ApplicationComponents/Resources/UI/pqExampleVisualizationsDialog.ui \
  < "${PATCH_DIR}/paraview-5.5.0-pqExampleVisualizationsDialog.ui.patch"
$PATCH_BIN Qt/ApplicationComponents/pqExampleVisualizationsDialog.cxx \
  < "${PATCH_DIR}/paraview-5.5.0-pqExampleVisualizationsDialog.cxx.patch"
$PATCH_BIN Remoting/Core/vtkPVFileInformation.cxx \
  < "${PATCH_DIR}/paraview-5.8.0-vtkPVFileInformation.cxx.patch"

## CPack variables for packaging meta-data
$PATCH_BIN -p1 \
  < "${PATCH_DIR}/paraview-5.9.0-CPack-CMakeLists.txt.patch"
$PATCH_BIN -p1 \
  < "${PATCH_DIR}/paraview-5.9.0-build-options-CMakeLists.txt.patch"
mkdir -p .github/workflows/
cp ${PATCH_DIR}/package.yml .github/workflows
cp ${PATCH_DIR}/headless.yml .github/workflows

echo "Finished patching."

cd "$PATCH_DIR" || exit 4
