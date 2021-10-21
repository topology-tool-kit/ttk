#include <iso646.h>
#include <regex>
#include <ttkImportEmbeddingFromTable.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// VTK includes
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkPointSet.h>
#include <vtkTable.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkImportEmbeddingFromTable)

  int ttkImportEmbeddingFromTable::FillInputPortInformation(
    int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }

  return 0;
}

int ttkImportEmbeddingFromTable::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename VTK_TT>
inline void setPointFromData(const vtkSmartPointer<vtkPoints> &points,
                             VTK_TT *xdata,
                             VTK_TT *ydata,
                             VTK_TT *zdata,
                             const bool Embedding2D) {
  for(SimplexId i = 0; i < points->GetNumberOfPoints(); ++i) {
    double p[3];
    p[0] = xdata[i];
    p[1] = ydata[i];
    p[2] = Embedding2D ? 0 : zdata[i];
    points->SetPoint(i, p);
  }
}

int ttkImportEmbeddingFromTable::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkPointSet *inputDataSet = vtkPointSet::GetData(inputVector[0]);
  vtkTable *inputTable = vtkTable::GetData(inputVector[1]);
  vtkPointSet *output = vtkPointSet::GetData(outputVector);

  const SimplexId numberOfPoints = inputDataSet->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPoints <= 0) {
    printErr("input has no point.");
    return -1;
  }
#endif

  vtkDataArray *xarr = XColumn.empty()
                         ? nullptr
                         : vtkDataArray::SafeDownCast(
                           inputTable->GetColumnByName(XColumn.data()));
  vtkDataArray *yarr = YColumn.empty()
                         ? nullptr
                         : vtkDataArray::SafeDownCast(
                           inputTable->GetColumnByName(YColumn.data()));
  vtkDataArray *zarr = ZColumn.empty()
                         ? nullptr
                         : vtkDataArray::SafeDownCast(
                           inputTable->GetColumnByName(ZColumn.data()));

#ifndef TTK_ENABLE_KAMIKAZE
  if(xarr == nullptr or yarr == nullptr or zarr == nullptr) {
    printErr("invalid input columns.");
    return -1;
  }
  if(xarr->GetNumberOfTuples() != numberOfPoints
     or yarr->GetNumberOfTuples() != numberOfPoints
     or zarr->GetNumberOfTuples() != numberOfPoints) {
    printErr("number of points on inputs mismatch.");
    return -1;
  }
  if(xarr->GetDataType() != yarr->GetDataType()
     or xarr->GetDataType() != zarr->GetDataType()) {
    printErr("input columns has different data types.");
    return -1;
  }
#endif

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numberOfPoints);

  switch(xarr->GetDataType()) {
    vtkTemplateMacro(setPointFromData(
      points, static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(xarr)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(yarr)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(zarr)), Embedding2D));
  }

  output->ShallowCopy(inputDataSet);
  output->SetPoints(points);

  return 1;
}
