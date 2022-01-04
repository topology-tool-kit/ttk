/// \ingroup vtk
/// \class ttkDepthImageBasedGeometryApproximation
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.07.2018
///
/// \brief TTK VTK-filter that approximates the geomerty that is depicted by a
/// set of depth images.
///
/// VTK wrapping code for the ttk::DepthImageBasedGeometryApproximation package.
///
/// This filter approximates the geometry that is depicted by a set of depth
/// images.
///
/// Related publication:
/// 'VOIDGA: A View-Approximation Oriented Image Database Generation Approach'
/// Jonas Lukasczyk, Eric Kinner, James Ahrens, Heike Leitte, and Christoph
/// Garth. IEEE 8th Symposium on Large Data Analysis and Visualization (LDAV),
/// 2018.
///
/// \param Input A vtkMultiBlockDataSet containing a set of depth images
/// represented by vtkImagedata objects. (vtkMultiBlockDataSet) \param
/// Subsampling The factor that controls the sampling rate (0: no Subsampling,
/// 1: skip every second sample, 2: skip every second and third sample...)
/// \param Output A set of unstructured grids where each grid corresponds to a
/// depth image (vtkMultiBlockDataSet)
///
/// \sa ttk::DepthImageBasedGeometryApproximation
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/geometryApproximation/">Geometry
///   Approximation example</a> \n

#pragma once

// VTK Module
#include <ttkDepthImageBasedGeometryApproximationModule.h>

// VTK includes
#include <ttkAlgorithm.h>

// TTK includes
#include <DepthImageBasedGeometryApproximation.h>

class TTKDEPTHIMAGEBASEDGEOMETRYAPPROXIMATION_EXPORT
  ttkDepthImageBasedGeometryApproximation
  : public ttkAlgorithm,
    protected ttk::DepthImageBasedGeometryApproximation {

public:
  static ttkDepthImageBasedGeometryApproximation *New();
  vtkTypeMacro(ttkDepthImageBasedGeometryApproximation, ttkAlgorithm);

protected:
  ttkDepthImageBasedGeometryApproximation();
  ~ttkDepthImageBasedGeometryApproximation() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
