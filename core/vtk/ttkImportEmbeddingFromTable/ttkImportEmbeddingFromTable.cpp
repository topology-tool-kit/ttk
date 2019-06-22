#include <iso646.h>
#include <regex>
#include <ttkImportEmbeddingFromTable.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkImportEmbeddingFromTable)

  // transmit abort signals
  bool ttkImportEmbeddingFromTable::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status
int ttkImportEmbeddingFromTable::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkImportEmbeddingFromTable] " << progress * 100
        << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkImportEmbeddingFromTable::doIt(vtkPointSet *inputDataSet,
                                      vtkTable *inputTable,
                                      vtkPointSet *output) {
  Memory m;

  const SimplexId numberOfPoints = inputDataSet->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPoints <= 0) {
    cerr << "[ttkImportEmbeddingFromTable] Error: input has no point." << endl;
    return -1;
  }
#endif

  vtkAbstractArray *xarr
    = XColumn.empty() ? nullptr : inputTable->GetColumnByName(XColumn.data());
  vtkAbstractArray *yarr
    = YColumn.empty() ? nullptr : inputTable->GetColumnByName(YColumn.data());
  vtkAbstractArray *zarr
    = ZColumn.empty() ? nullptr : inputTable->GetColumnByName(ZColumn.data());

#ifndef TTK_ENABLE_KAMIKAZE
  if(xarr == nullptr or yarr == nullptr or zarr == nullptr) {
    cerr << "[ttkImportEmbeddingFromTable] Error: invalid input columns."
         << endl;
    return -1;
  }
  if(xarr->GetNumberOfTuples() != numberOfPoints
     or yarr->GetNumberOfTuples() != numberOfPoints
     or zarr->GetNumberOfTuples() != numberOfPoints) {
    cerr << "[ttkImportEmbeddingFromTable] Error: number of points on inputs "
            "mismatch."
         << endl;
    return -1;
  }
  if(xarr->GetDataType() != yarr->GetDataType()
     or xarr->GetDataType() != zarr->GetDataType()) {
    cerr << "[ttkImportEmbeddingFromTable] Error: input columns has different "
            "data types."
         << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numberOfPoints);

  switch(xarr->GetDataType()) {
    vtkTemplateMacro({
      VTK_TT *xdata = static_cast<VTK_TT *>(xarr->GetVoidPointer(0));
      VTK_TT *ydata = static_cast<VTK_TT *>(yarr->GetVoidPointer(0));
      VTK_TT *zdata = static_cast<VTK_TT *>(zarr->GetVoidPointer(0));

      for(SimplexId i = 0; i < numberOfPoints; ++i) {
        double p[3];
        p[0] = xdata[i];
        p[1] = ydata[i];
        p[2] = Embedding2D ? 0 : zdata[i];
        points->SetPoint(i, p);
      }
    });
  }

  output->ShallowCopy(inputDataSet);
  output->SetPoints(points);

  {
    stringstream msg;
    msg << "[ttkImportEmbeddingFromTable] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkImportEmbeddingFromTable::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  Memory m;

  vtkPointSet *inputDataSet = vtkPointSet::GetData(inputVector[0]);
  vtkTable *inputTable = vtkTable::GetData(inputVector[1]);
  vtkPointSet *output = vtkPointSet::GetData(outputVector);

  doIt(inputDataSet, inputTable, output);

  {
    stringstream msg;
    msg << "[ttkImportEmbeddingFromTable] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
