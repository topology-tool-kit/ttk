#include <ttkDataSetToTable.h>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>

vtkStandardNewMacro(ttkDataSetToTable);

ttkDataSetToTable::ttkDataSetToTable() {
  this->setDebugMsgPrefix("DataSetToTable");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkDataSetToTable::~ttkDataSetToTable() {
}

int ttkDataSetToTable::FillInputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkDataSetToTable::FillOutputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDataSetToTable::RequestData(vtkInformation *ttkNotUsed(request),
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  ttk::Timer t;

  std::string targetAttributeName = this->DataAssociation == 0   ? "Point"
                                    : this->DataAssociation == 1 ? "Cell"
                                                                 : "Field";

  this->printMsg("Converting " + targetAttributeName + "Data to Table", 0, 0,
                 ttk::debug::LineMode::REPLACE);

  auto input = vtkDataSet::GetData(inputVector[0]);
  if(!input)
    return 0;

  auto targetData = this->DataAssociation == 0   ? input->GetPointData()
                    : this->DataAssociation == 1 ? input->GetCellData()
                                                 : input->GetFieldData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(targetData->GetNumberOfArrays() < 1) {
    this->printWrn("Empty " + targetAttributeName + "Data.");
  }
#endif

  auto table = vtkTable::GetData(outputVector);
  table->GetRowData()->ShallowCopy(targetData);

  this->printMsg("Converting " + targetAttributeName + "Data to Table", 1,
                 t.getElapsedTime());

  return 1;
}