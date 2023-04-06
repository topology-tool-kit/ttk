#include <ttkTimeVaryingPersistenceDiagramDistanceMatrix.h>
#include <ttkPersistenceDiagramUtils.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTimeVaryingPersistenceDiagramDistanceMatrix);

ttkTimeVaryingPersistenceDiagramDistanceMatrix::ttkTimeVaryingPersistenceDiagramDistanceMatrix() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkTimeVaryingPersistenceDiagramDistanceMatrix::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkTimeVaryingPersistenceDiagramDistanceMatrix::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0; 
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkTimeVaryingPersistenceDiagramDistanceMatrix::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Memory m;

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  int nBlocks = blocks->GetNumberOfBlocks();
  cout << "nBlocks in the first vtkMultiBlockDataSet (inputVector[0]): " << nBlocks << endl ;
    
    for(int i = 0; i < nBlocks; i++){
      auto block = blocks->GetBlock(i);
      vtkMultiBlockDataSet * multiBlockDataSet = vtkMultiBlockDataSet::SafeDownCast(block);
      cout << "block " << i <<" has "<< multiBlockDataSet->GetNumberOfBlocks() << " elements" << endl;
  }
  
  return 1;
}
