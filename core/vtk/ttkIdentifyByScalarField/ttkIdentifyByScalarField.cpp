#include <numeric>
#include <ttkIdentifyByScalarField.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifyByScalarField)

  int ttkIdentifyByScalarField::getScalars(vtkDataSet *input) {
  vtkCellData *cellData = input->GetCellData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellData) {
    cerr << "[ttkIdentifyByScalarField] Error : input has no point data."
         << endl;
    return -1;
  }

  if(!ScalarField.length()) {
    cerr << "[ttkIdentifyByScalarField] Error : scalar field has no name."
         << endl;
    return -2;
  }
#endif

  if(ScalarField.length()) {
    inputScalars_ = cellData->GetArray(ScalarField.data());
  } else {
    inputScalars_ = cellData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_) {
    cerr << "[ttkIdentifyByScalarField] Error : input scalar field pointer is "
            "null."
         << endl;
    return -3;
  }
#endif

  return 0;
}

int ttkIdentifyByScalarField::doIt(vector<vtkDataSet *> &inputs,
                                   vector<vtkDataSet *> &outputs) {
  Memory m;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs.size()) {
    cerr << "[ttkIdentifyByScalarField] Error: not enough input information."
         << endl;
    return -1;
  }
#endif

  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkIdentifyByScalarField] Error: input pointer is NULL." << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkIdentifyByScalarField] Error: input has no point." << endl;
    return -1;
  }

  if(!input->GetNumberOfCells()) {
    cerr << "[ttkIdentifyByScalarField] Error: input has no cell." << endl;
    return -1;
  }
#endif

  if(getScalars(input)) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkIndentifyByScalarField] Error : wrong scalars." << endl;
    return -1;
#endif
  }

  const SimplexId numberOfCells = input->GetNumberOfCells();
  vector<SimplexId> inputIds(numberOfCells);
  std::iota(inputIds.begin(), inputIds.end(), 0);
  switch(inputScalars_->GetDataType()) {
    ttkTemplateMacro({
      VTK_TT *arr = static_cast<VTK_TT *>(inputScalars_->GetVoidPointer(0));

      auto greater_cmp = [arr](int a, int b) { return arr[a] > arr[b]; };

      auto lower_cmp = [arr](int a, int b) { return arr[a] < arr[b]; };

      if(IncreasingOrder)
        std::sort(inputIds.begin(), inputIds.end(), lower_cmp);
      else
        std::sort(inputIds.begin(), inputIds.end(), greater_cmp);
    });
  }

  vtkSmartPointer<ttkSimplexIdTypeArray> ids
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  ids->SetNumberOfComponents(1);
  ids->SetNumberOfTuples(numberOfCells);
  ids->SetName("CellScalarFieldName");

  SimplexId *outputIds = static_cast<SimplexId *>(ids->GetVoidPointer(0));

  for(int i = 0; i < numberOfCells; ++i)
    outputIds[inputIds[i]] = i;
  if(StartByOne) {
    for(int i = 0; i < numberOfCells; ++i)
      outputIds[i] += 1;
  }

  output->ShallowCopy(input);

  output->GetCellData()->AddArray(ids);

  {
    stringstream msg;
    msg << "[ttkIdentifyByScalarField] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
