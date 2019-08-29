#include <ttkDataSetToTable.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDataSetToTable)

  // transmit abort signals
  bool ttkDataSetToTable::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status
int ttkDataSetToTable::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkDataSetToTable] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkDataSetToTable::doIt(vtkDataSet *input, vtkTable *output) {
  Memory m;

  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

  switch(DataAssociation) {
    case AssociationType::Point: {
      vtkPointData *inputPointData = input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!inputPointData) {
        cerr << "[ttkDataSetToTable] Error: input has no point data." << endl;
        return -1;
      }
#endif

      const SimplexId numberOfArrays = inputPointData->GetNumberOfArrays();
#ifndef TTK_ENABLE_KAMIKAZE
      if(numberOfArrays <= 0) {
        cerr << "[ttkDataSetToTable] Error: input point data is empty." << endl;
        return -1;
      }
#endif

      table->GetRowData()->ShallowCopy(inputPointData);
    } break;

    case AssociationType::Cell: {
      vtkCellData *inputCellData = input->GetCellData();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!inputCellData) {
        cerr << "[ttkDataSetToTable] Error: input has no cell data." << endl;
        return -1;
      }
#endif

      const SimplexId numberOfArrays = inputCellData->GetNumberOfArrays();
#ifndef TTK_ENABLE_KAMIKAZE
      if(numberOfArrays <= 0) {
        cerr << "[ttkDataSetToTable] Error: input cell data is empty." << endl;
        return -1;
      }
#endif

      table->GetRowData()->ShallowCopy(inputCellData);
    } break;
  }

  output->ShallowCopy(table);

  {
    stringstream msg;
    msg << "[ttkDataSetToTable] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkDataSetToTable::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  Memory m;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkTable *output = vtkTable::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkDataSetToTable] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
