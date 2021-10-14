/// \ingroup vtk
/// \class ttkTriangulationRequest
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date April 2017.
///
/// \brief TTK VTK-filter that wraps the triangulation processing package.
///
/// VTK wrapping code for the @Triangulation package.
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

// VTK Module
#include <ttkTriangulationRequestModule.h>

class TTKTRIANGULATIONREQUEST_EXPORT ttkTriangulationRequest
  : public ttkAlgorithm {

public:
  enum Simplex { Vertex = 0, Edge, Triangle, Tetra };

  enum Request {
    ComputeSimplex = 0,
    ComputeFacet,
    ComputeCofacet,
    ComputeStar,
    ComputeLink
  };

  static ttkTriangulationRequest *New();
  vtkTypeMacro(ttkTriangulationRequest, ttkAlgorithm)

    vtkSetMacro(SimplexType, int);
  vtkGetMacro(SimplexType, int);

  vtkSetMacro(SimplexIdentifier, int);
  vtkGetMacro(SimplexIdentifier, int);

  vtkSetMacro(RequestType, int);
  vtkGetMacro(RequestType, int);

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

protected:
  ttkTriangulationRequest() {
    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
    this->setDebugMsgPrefix("TriangulationRequest");
  }

  ~ttkTriangulationRequest() override{};

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int SimplexType{0};
  int SimplexIdentifier{0};
  int RequestType{0};
  bool KeepAllDataArrays{true};
};
