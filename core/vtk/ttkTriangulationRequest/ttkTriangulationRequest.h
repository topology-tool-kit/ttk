/// \ingroup vtk
/// \class ttkTriangulationRequest
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date April 2017.
///
/// \brief TTK VTK-filter that wraps the triangulation processing package.
///
/// VTK wrapping code for the ttk::Triangulation package.
///
/// \param Input Geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output Output geometry requested by the user (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::Triangulation
///
#pragma once

// ttk code includes
#include <Triangulation.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

// VTK Module
#include <ttkTriangulationRequestModule.h>

class TTKTRIANGULATIONREQUEST_EXPORT ttkTriangulationRequest
  : public ttkAlgorithm {

public:
  enum class SIMPLEX {
    VERTEX = 0,
    EDGE = 1,
    TRIANGLE = 2,
    TETRA = 3,
  };
  enum class REQUEST {
    COMPUTE_SIMPLEX = 0,
    COMPUTE_FACET = 1,
    COMPUTE_COFACET = 2,
    COMPUTE_STAR = 3,
    COMPUTE_LINK = 4,
  };

  static ttkTriangulationRequest *New();
  vtkTypeMacro(ttkTriangulationRequest, ttkAlgorithm);

  ttkSetEnumMacro(SimplexType, SIMPLEX);
  vtkGetEnumMacro(SimplexType, SIMPLEX);

  vtkSetMacro(SimplexIdentifier, const std::string &);
  vtkGetMacro(SimplexIdentifier, std::string);

  ttkSetEnumMacro(RequestType, REQUEST);
  vtkGetEnumMacro(RequestType, REQUEST);

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

protected:
  ttkTriangulationRequest();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  SIMPLEX SimplexType{};
  REQUEST RequestType{};
  std::string SimplexIdentifier{0};
  bool KeepAllDataArrays{true};
};
