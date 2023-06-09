/// \ingroup vtk
/// \class ttkMarchingTetrahedra
/// \author Robin Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief TTK VTK-filter that wraps the ttk::MarchingTetrahedra module.
///
/// Given an input point data array and triangulation this class executes the
/// marching tetrahedra/triangles algorithm. It has three options that either
/// separate each label with a single separating geometry inbetween two labels,
/// or a separating geometry enclosing each label (detailed and fast mode).
///
/// \param Input Input scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output Output separating geometry (vtkPolyData)
///
/// This filter can be used like any other VTK filter (for instance, by using
/// the sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MarchingTetrahedra/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Parallel Computation of Piecewise Linear Morse-Smale Segmentations" \n
/// Robin G. C. Maack, Jonas Lukasczyk, Julien Tierny, Hans Hagen,
/// Ross Maciejewski, Christoph Garth \n
/// IEEE Transactions on Visualization and Computer Graphics \n
///
/// \sa ttk::MarchingTetrahedra
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleSegmentation_at/">Morse-Smale
///   segmentation example</a> \n

#pragma once

// VTK Module
#include <ttkMarchingTetrahedraModule.h>

// ttk code includes
#include <MarchingTetrahedra.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class vtkPolyData;

class TTKMARCHINGTETRAHEDRA_EXPORT ttkMarchingTetrahedra
  : public ttkAlgorithm,
    protected ttk::MarchingTetrahedra {

public:
  static ttkMarchingTetrahedra *New();

  vtkTypeMacro(ttkMarchingTetrahedra, ttkAlgorithm);

  ttkSetEnumMacro(SurfaceMode, SURFACE_MODE);
  vtkGetEnumMacro(SurfaceMode, SURFACE_MODE);

protected:
  template <typename scalarType, typename triangulationType>
  int dispatch(vtkDataArray *const inputScalars,
               vtkPolyData *const outputSeparators,
               const triangulationType &triangulation);

  ttkMarchingTetrahedra();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
