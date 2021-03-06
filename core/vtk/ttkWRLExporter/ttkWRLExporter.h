/// \ingroup vtk
/// \class ttkWRLExporter
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK helper that fixes a few bugs in the texture support for the VRML
/// export.

#pragma once

#include <ttkWRLExporterModule.h>

class vtkPolyData;
TTKWRLEXPORTER_EXPORT vtkPolyData *wrlExporterPolyData_;
