#include <ttkIdentifiers.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <map>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

vtkCellArray *GetCells(vtkDataSet *dataSet) {
  switch(dataSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID: {
      auto dataSetAsUG = static_cast<vtkUnstructuredGrid *>(dataSet);
      return dataSetAsUG->GetCells();
    }
    case VTK_POLY_DATA: {
      auto dataSetAsPD = static_cast<vtkPolyData *>(dataSet);
      return dataSetAsPD->GetNumberOfPolys() > 0   ? dataSetAsPD->GetPolys()
             : dataSetAsPD->GetNumberOfLines() > 0 ? dataSetAsPD->GetLines()
                                                   : dataSetAsPD->GetVerts();
    }
  }
  return nullptr;
}

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifiers);

ttkIdentifiers::ttkIdentifiers() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  setDebugMsgPrefix("Identifiers");

  //   vtkWarningMacro("`TTK Identifiers' is now deprecated. Please use "
  //                   "`Generate Ids' instead.");
}

ttkIdentifiers::~ttkIdentifiers() = default;

int ttkIdentifiers::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkIdentifiers::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkIdentifiers::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);

  Timer t;

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Triangulation is NULL");
    return 0;
  }

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
