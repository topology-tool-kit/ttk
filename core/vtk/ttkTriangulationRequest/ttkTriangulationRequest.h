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
#include <ttkWrapper.h>

// VTK includes
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTriangulationRequest
#else
class ttkTriangulationRequest
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

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
  vtkTypeMacro(ttkTriangulationRequest, vtkDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }
  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  vtkSetMacro(SimplexType, int);
  vtkGetMacro(SimplexType, int);

  vtkSetMacro(SimplexIdentifier, int);
  vtkGetMacro(SimplexIdentifier, int);

  vtkSetMacro(RequestType, int);
  vtkGetMacro(RequestType, int);

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

  int FillInputPortInformation(int port, vtkInformation *info) override {

    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
        break;

      default:
        break;
    }

    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {

    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;

      default:
        break;
    }

    return 1;
  }

protected:
  ttkTriangulationRequest() {
    UseAllCores = true;
    SimplexType = 0;
    SimplexIdentifier = 0;
    RequestType = 0;
    KeepAllDataArrays = true;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  ~ttkTriangulationRequest(){};

  TTK_SETUP();

private:
  int SimplexType;
  int SimplexIdentifier;
  int RequestType;
  bool KeepAllDataArrays;
};
