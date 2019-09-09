/// \ingroup vtk
/// \class ttkMeshGraph
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.12.2018
///
/// \brief TTK VTK-filter that generates a mesh for a graph.
///
/// This filter generates for each one dimensional cell (edge) of a
/// 'vtkUnstructuredGrid' a two dimensional cell by mapping a size value to the
/// width of the input cell. The output is a 'vtkUnstructuredGrid' consisting of
/// a set of either quadratic quads or linear polygons.
///
/// VTK wrapping code for the @MeshGraph package.
///
/// \param Input Graph (vtkUnstructuredGrid)
/// \param Output Graph (vtkUnstructuredGrid)
///
/// \sa ttk::MeshGraph

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkUnstructuredGridAlgorithm.h>

// TTK includes
#include <MeshGraph.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkMeshGraph
#else
class ttkMeshGraph
#endif
  : public vtkUnstructuredGridAlgorithm,
    public ttk::Wrapper {

public:
  static ttkMeshGraph *New();
  vtkTypeMacro(ttkMeshGraph, vtkUnstructuredGridAlgorithm)

    vtkSetMacro(UseVariableSize, bool);
  vtkGetMacro(UseVariableSize, bool);

  vtkSetMacro(SizeFieldName, std::string);
  vtkGetMacro(SizeFieldName, std::string);

  vtkSetMacro(SizeAxis, int);
  vtkGetMacro(SizeAxis, int);

  vtkSetMacro(SizeScale, float);
  vtkGetMacro(SizeScale, float);

  vtkSetMacro(UseQuadraticCells, bool);
  vtkGetMacro(UseQuadraticCells, bool);

  vtkSetMacro(Subdivisions, int);
  vtkGetMacro(Subdivisions, int);

  vtkSetMacro(Tetrahedralize, bool);
  vtkGetMacro(Tetrahedralize, bool);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);
  void SetThreads() {
    threadNumber_
      = !UseAllCores ? ThreadNumber : ttk::OsCall::getNumberOfCores();
    Modified();
  }
  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }
  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      default:
        return 0;
    }
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkMeshGraph() {

    SetUseVariableSize(true);
    SetSizeFieldName("Size");

    SetSizeAxis(0);
    SetSizeScale(1);
    SetUseQuadraticCells(true);
    SetSubdivisions(0);
    SetTetrahedralize(true);

    UseAllCores = false;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkMeshGraph(){};

  bool UseAllCores;
  int ThreadNumber;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseVariableSize;
  std::string SizeFieldName;

  int SizeAxis;
  float SizeScale;
  bool UseQuadraticCells;
  int Subdivisions;
  bool Tetrahedralize;

  ttk::MeshGraph meshGraph;

  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};