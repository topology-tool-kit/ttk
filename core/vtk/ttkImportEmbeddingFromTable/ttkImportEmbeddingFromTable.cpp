#include <ttkImportEmbeddingFromTable.h>

#include <ttkUtils.h>
#include <vtkPointSet.h>
#include <vtkTable.h>

#include <iso646.h>
#include <regex>

vtkStandardNewMacro(ttkImportEmbeddingFromTable);

ttkImportEmbeddingFromTable::ttkImportEmbeddingFromTable() {
  this->setDebugMsgPrefix("EmbeddingFromTable");

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkImportEmbeddingFromTable::~ttkImportEmbeddingFromTable() {
}

int ttkImportEmbeddingFromTable::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if (port == 1){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkImportEmbeddingFromTable::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename VTK_TT>
inline void setPointFromData(vtkSmartPointer<vtkPoints> points,
                             VTK_TT *xdata,
                             VTK_TT *ydata,
                             VTK_TT *zdata,
                             const bool Embedding2D) {
  for(size_t i = 0; i < points->GetNumberOfPoints(); ++i) {
    double p[3];
    p[0] = xdata[i];
    p[1] = ydata[i];
    p[2] = Embedding2D ? 0 : zdata[i];
    points->SetPoint(i, p);
  }
}

int ttkImportEmbeddingFromTable::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  ttk::Timer t;

  this->printMsg(
    "Apply Embedding",
    0, 0,
    ttk::debug::LineMode::REPLACE
  );

  auto inputDataSet = vtkPointSet::GetData(inputVector[0]);
  auto inputTable = vtkTable::GetData(inputVector[1]);
  auto outputDataSet = vtkPointSet::GetData(outputVector);

  const size_t numberOfPoints = inputDataSet->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPoints <= 0) {
    this->printErr("Input has no points.");
    return 0;
  }
#endif

  auto xarr = this->GetInputArrayToProcess(0, inputVector);
  auto yarr = this->GetInputArrayToProcess(1, inputVector);
  auto zarr = this->GetInputArrayToProcess(2, inputVector);

#ifndef TTK_ENABLE_KAMIKAZE
  if(xarr == nullptr or yarr == nullptr or zarr == nullptr) {
    this->printErr("Invalid input columns.");
    return 0;
  }
  if(xarr->GetNumberOfTuples() != numberOfPoints
     or yarr->GetNumberOfTuples() != numberOfPoints
     or zarr->GetNumberOfTuples() != numberOfPoints) {
    this->printErr("Number of points on inputs mismatch.");
    return 0;
  }
  if(xarr->GetDataType() != yarr->GetDataType()
     or xarr->GetDataType() != zarr->GetDataType()) {
    this->printErr("Input columns has different data types.");
    return 0;
  }
#endif

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numberOfPoints);

  switch(xarr->GetDataType()) {
    vtkTemplateMacro(
      setPointFromData(
        points,
        (VTK_TT *) ttkUtils::GetVoidPointer(xarr),
        (VTK_TT *) ttkUtils::GetVoidPointer(yarr),
        (VTK_TT *) ttkUtils::GetVoidPointer(zarr),
        Embedding2D
      )
    );
  }

  outputDataSet->ShallowCopy(inputDataSet);
  outputDataSet->SetPoints(points);

  this->printMsg(
    "Apply Embedding",
    1, t.getElapsedTime()
  );

  return 1;
}