/// \ingroup vtk
/// \class ttkLDistanceMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief Computes a distance matrix using LDistance between several
/// same-dimensions input datasets
///
/// \sa LDistanceMatrix

#pragma once

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkLDistanceMatrixModule.h>
//
#include <ttkTriangulationAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKLDISTANCEMATRIX_EXPORT ttkLDistanceMatrix
  : public vtkMultiBlockDataSetAlgorithm,
    protected ttk::Wrapper {

public:
  static ttkLDistanceMatrix *New();

  vtkTypeMacro(ttkLDistanceMatrix, vtkMultiBlockDataSetAlgorithm);

  // default ttk setters
  void SetDebugLevel(int debugLevel) {
    setDebugLevel(debugLevel);
    Modified();
  }

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  void SetThreadNumber(int data) {
    ThreadNumber = data;
    Modified();
  }
  // end of default ttk setters

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(DistanceType, std::string);
  vtkGetMacro(DistanceType, std::string);

protected:
  ttkLDistanceMatrix();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  bool needsToAbort() override {
    return GetAbortExecute();
  }
  int updateProgress(const float &progress) override {
    return 0;
  }

private:
  std::string DistanceType{};
  std::string ScalarField{};
  int ThreadNumber{1};
  bool UseAllCores{true};
};
