#include <ttkCinemaImaging.h>

#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

#include <ttkUtils.h>
#include <vtkDoubleArray.h>

#include <ttkCinemaImagingEmbree.h>
#include <ttkCinemaImagingNative.h>
#include <ttkCinemaImagingVTK.h>

#include <CinemaImaging.h>

vtkStandardNewMacro(ttkCinemaImaging);

ttkCinemaImaging::ttkCinemaImaging() {
  this->setDebugMsgPrefix("CinemaImaging");
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
};

ttkCinemaImaging::~ttkCinemaImaging(){};

int ttkCinemaImaging::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  } else if(port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;
  return 1;
};

int ttkCinemaImaging::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
};

int ttkCinemaImaging::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  auto inputObject = vtkDataObject::GetData(inputVector[0]);
  auto inputGrid = vtkPointSet::GetData(inputVector[1]);
  auto outputCollection = vtkMultiBlockDataSet::GetData(outputVector);
  if(!inputObject || !inputGrid) {
    this->printErr("Unsupported input object types.");
    return 0;
  }

  auto inputObjectAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(inputObject->IsA("vtkMultiBlockDataSet"))
    inputObjectAsMB->ShallowCopy(inputObject);
  else
    inputObjectAsMB->SetBlock(0, inputObject);
  const size_t nInputObjects = inputObjectAsMB->GetNumberOfBlocks();

  // get default parameters
  std::vector<double> defaultFocalPoint{
    this->FocalPoint[0], this->FocalPoint[1], this->FocalPoint[2]};
  std::vector<double> defaultNearFar{this->NearFar[0], this->NearFar[1]};
  double defaultHeight = this->Height;
  double defaultAngle = this->Angle;

  if(this->AutoFocalPoint || this->AutoNearFar || this->AutoHeight) {

    double objectBounds[6];
    inputObjectAsMB->GetBounds(objectBounds);

    auto norm = [](const double v[3]) {
      return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };

    const double d[3]{objectBounds[1] - objectBounds[0],
                      objectBounds[3] - objectBounds[2],
                      objectBounds[5] - objectBounds[4]};
    const double objectDiameter = norm(d);

    const double c[3]{objectBounds[0] + 0.5 * d[0],
                      objectBounds[2] + 0.5 * d[1],
                      objectBounds[4] + 0.5 * d[2]};

    if(this->AutoFocalPoint) {
      defaultFocalPoint[0] = c[0];
      defaultFocalPoint[1] = c[1];
      defaultFocalPoint[2] = c[2];
    }

    if(this->AutoNearFar) {
      double gridBounds[6];
      inputGrid->GetBounds(gridBounds);

      const double l[3]{
        gridBounds[0] - c[0],
        gridBounds[2] - c[1],
        gridBounds[4] - c[2],
      };
      const double ld = norm(l);

      defaultNearFar[0] = std::max(0.01, ld - objectDiameter * 0.5);
      defaultNearFar[1] = ld + objectDiameter * 0.5;
    }

    if(this->AutoHeight) {
      defaultHeight = objectDiameter;
    }
  }

  // ensure grid contains all camera parameters as field data
  auto aInputGrid
    = vtkSmartPointer<vtkPointSet>::Take(inputGrid->NewInstance());
  aInputGrid->ShallowCopy(inputGrid);
  int n = aInputGrid->GetNumberOfPoints();
  auto aInputGridPD = aInputGrid->GetPointData();

  ttkCinemaImaging::EnsureGridData(
    aInputGridPD, "Resolution", n,
    {(double)this->Resolution[0], (double)this->Resolution[1]});
  ttkCinemaImaging::EnsureGridData(
    aInputGridPD, "ProjectionMode", n, {(double)this->ProjectionMode});

  ttkCinemaImaging::EnsureGridData(
    aInputGridPD, "CamNearFar", n, defaultNearFar);
  ttkCinemaImaging::EnsureGridData(
    aInputGridPD, "CamHeight", n, {defaultHeight});
  ttkCinemaImaging::EnsureGridData(aInputGridPD, "CamAngle", n, {defaultAngle});
  ttkCinemaImaging::EnsureGridData(aInputGridPD, "CamUp", n, {0, 1, 0});

  const bool gridHasDirection = aInputGridPD->HasArray("CamDirection");
  const bool gridHasFocalPoint = aInputGridPD->HasArray("CamFocalPoint");
  if(gridHasFocalPoint || !gridHasDirection) {
    ttkCinemaImaging::EnsureGridData(
      aInputGridPD, "CamFocalPoint", n, defaultFocalPoint);
    ttkCinemaImaging::ComputeDirFromFocalPoint(aInputGrid);
  } else {
    ttkCinemaImaging::EnsureGridData(
      aInputGridPD, "CamDirection", n, {0, 0, -1});
  }

  // print parameters
  {
    this->printMsg(
      {{"#Objects", std::to_string(nInputObjects)},
       {"Backend", this->Backend == 0   ? std::string("VTK_OPENGL")
                   : this->Backend == 1 ? std::string("EMBREE")
                                        : std::string("NATIVE")},
       {"Resolution", std::to_string(this->Resolution[0]) + " x "
                        + std::to_string(this->Resolution[1])},
       {"Projection", this->ProjectionMode ? "Perspective" : "Orthographic"},
       {"CamNearFar", std::to_string(defaultNearFar[0]) + " / "
                        + std::to_string(defaultNearFar[1])},
       {"CamFocalPoint", gridHasDirection || gridHasFocalPoint
                           ? "Grid Data"
                           : "(" + std::to_string(defaultFocalPoint[0]) + ", "
                               + std::to_string(defaultFocalPoint[1]) + ", "
                               + std::to_string(defaultFocalPoint[2]) + ")"},
       {std::string(this->ProjectionMode == 0 ? "CamHeight" : "CamAngle"),
        std::to_string(this->ProjectionMode == 0 ? defaultHeight
                                                 : defaultAngle)}});
    this->printMsg(ttk::debug::Separator::L1);
  }

  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  for(size_t b = 0; b < nInputObjects; b++) {
    auto outputImages = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    auto inputPointSet
      = vtkPointSet::SafeDownCast(inputObjectAsMB->GetBlock(b));
    if(!inputPointSet->IsA("vtkUnstructuredGrid")
       && !inputPointSet->IsA("vtkPolyData")) {
      this->printErr("Unsupported input object type.");
      return 0;
    }
    this->RequestDataSingle(outputImages,

                            inputPointSet, aInputGrid, defaultFocalPoint,
                            defaultNearFar, defaultHeight, defaultAngle);
    outputAsMB->SetBlock(b, outputImages);
  }

  if(inputObject->IsA("vtkMultiBlockDataSet"))
    outputCollection->ShallowCopy(outputAsMB);
  else
    outputCollection->ShallowCopy(outputAsMB->GetBlock(0));

  return 1;
}

int ttkCinemaImaging::RequestDataSingle(
  vtkMultiBlockDataSet *outputImages,

  vtkPointSet *inputObject,
  vtkPointSet *inputGrid,
  const std::vector<double> &ttkNotUsed(defaultFocalPoint),
  const std::vector<double> &ttkNotUsed(defaultNearFar),
  const double ttkNotUsed(defaultHeight),
  const double ttkNotUsed(defaultAngle)) {
  ttk::Timer globalTimer;

  auto cells = ttkCinemaImaging::GetCells(inputObject);
  if(!cells)
    return 0;

  size_t nTriangles = cells->GetNumberOfCells();
  // make sure that cells consists only of triangles
  {
    auto offsets = static_cast<vtkIdType *>(
      ttkUtils::GetVoidPointer(cells->GetOffsetsArray()));
    vtkIdType j = 0;
    for(size_t i = 0; i < nTriangles; i++, j += 3) {
      if(j != offsets[i]) {
        this->printErr("Connectivity list has to contain only triangles.");
        return 0;
      }
    }
  }

  int status = 0;
  if(this->Backend == 0) {
    ttk::ttkCinemaImagingVTK renderer;
    renderer.setDebugLevel(this->debugLevel_);
    renderer.setThreadNumber(this->threadNumber_);
    status = renderer.RenderVTKObject(outputImages, inputObject, inputGrid);
  } else if(this->Backend == 1) {
    ttk::ttkCinemaImagingEmbree renderer;
    renderer.setDebugLevel(this->debugLevel_);
    renderer.setThreadNumber(this->threadNumber_);
    status = renderer.RenderVTKObject(outputImages, inputObject, inputGrid);
  } else {
    ttk::ttkCinemaImagingNative renderer;
    renderer.setDebugLevel(this->debugLevel_);
    renderer.setThreadNumber(this->threadNumber_);
    status = renderer.RenderVTKObject(outputImages, inputObject, inputGrid);
  }

  if(!status)
    return 0;

  return 1;
}

vtkCellArray *ttkCinemaImaging::GetCells(vtkPointSet *pointSet) {
  switch(pointSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID:
      return static_cast<vtkUnstructuredGrid *>(pointSet)->GetCells();
    case VTK_POLY_DATA:
      return static_cast<vtkPolyData *>(pointSet)->GetPolys();
  }
  return nullptr;
};

int ttkCinemaImaging::AddFieldDataArray(vtkFieldData *fd,
                                        vtkDataArray *array,
                                        int tupelIdx,
                                        const std::string &name) {
  if(!array)
    return 0;

  size_t nComponents = array->GetNumberOfComponents();

  auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
  newArray->SetName(name.empty() ? array->GetName() : name.data());
  newArray->SetNumberOfComponents(nComponents);
  newArray->SetNumberOfTuples(1);

  if(newArray->GetDataType() == array->GetDataType()) {
    newArray->SetTuple(0, tupelIdx, array);
  } else {
    for(size_t i = 0; i < nComponents; i++)
      newArray->SetValue(
        i, array->GetVariantValue(tupelIdx * nComponents + i).ToDouble());
  }

  fd->AddArray(newArray);

  return 1;
};

int ttkCinemaImaging::AddAllFieldDataArrays(vtkPointSet *inputGrid,
                                            vtkImageData *image,
                                            int tupelIdx) {
  auto imageFD = image->GetFieldData();

  auto inputGridPD = inputGrid->GetPointData();
  for(int i = 0; i < inputGridPD->GetNumberOfArrays(); i++) {
    ttkCinemaImaging::AddFieldDataArray(
      imageFD, inputGridPD->GetArray(i), tupelIdx);
  }

  ttkCinemaImaging::AddFieldDataArray(
    imageFD, inputGrid->GetPoints()->GetData(), tupelIdx, "CamPosition");

  return 1;
};

int ttkCinemaImaging::ComputeDirFromFocalPoint(vtkPointSet *inputGrid) {
  auto pos
    = static_cast<float *>(ttkUtils::GetVoidPointer(inputGrid->GetPoints()));
  auto focal = static_cast<double *>(ttkUtils::GetVoidPointer(
    inputGrid->GetPointData()->GetArray("CamFocalPoint")));

  int nTuples = inputGrid->GetNumberOfPoints();

  auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
  newArray->SetName("CamDirection");
  newArray->SetNumberOfComponents(3);
  newArray->SetNumberOfTuples(nTuples);

  auto dir = static_cast<double *>(ttkUtils::GetVoidPointer(newArray));

  for(int i = 0, j = nTuples * 3; i < j; i++)
    dir[i] = focal[i] - pos[i];

  inputGrid->GetPointData()->AddArray(newArray);

  return 1;
};

int ttkCinemaImaging::EnsureGridData(vtkPointData *fd,
                                     const std::string &name,
                                     int nTuples,
                                     const std::vector<double> &defaultValues) {
  auto array = vtkDoubleArray::SafeDownCast(fd->GetArray(name.data()));

  if(!array) {
    int nComponents = defaultValues.size();

    auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
    newArray->SetName(name.data());
    newArray->SetNumberOfComponents(nComponents);
    newArray->SetNumberOfTuples(nTuples);

    auto data = static_cast<double *>(ttkUtils::GetVoidPointer(newArray));
    for(int i = 0, q = 0; i < nTuples; i++)
      for(int j = 0; j < nComponents; j++)
        data[q++] = defaultValues[j];

    fd->AddArray(newArray);
  }

  return 1;
};

int ttkCinemaImaging::Normalize(vtkDataArray *depthArray,
                                const double nearFar[2]) {
  if(!depthArray->IsA("vtkFloatArray")
     || depthArray->GetNumberOfComponents() != 1)
    return 0;

  if(nearFar[0] == 0.0 && nearFar[1] == 0.0)
    return 1;

  const size_t nValues = depthArray->GetNumberOfTuples();
  auto data = static_cast<float *>(ttkUtils::GetVoidPointer(depthArray));

  const float near = (float)nearFar[0];
  const float far = (float)nearFar[1];
  const float delta = far - near;

  for(size_t i = 0; i < nValues; i++) {
    if(std::isnan(data[i])) {
      data[i] = 1.0;
    } else {
      data[i] = (data[i] - near) / delta;
      if(data[i] > 1.0)
        data[i] = 1.0;
      else if(data[i] < 0.0)
        data[i] = 0.0;
    }
  }

  return 1;
};

int ttkCinemaImaging::MapPointAndCellData(
  vtkImageData *outputImage,

  vtkPointSet *inputObject,
  const ttk::CinemaImaging *renderer,
  const unsigned int *primitiveIdArray,
  const float *barycentricCoordinates,
  const vtkIdType *inputObjectConnectivityList) {

  auto inputObjectPD = inputObject->GetPointData();
  auto inputObjectCD = inputObject->GetCellData();
  auto outputImagePD = outputImage->GetPointData();
  int dim[3];
  outputImage->GetDimensions(dim);
  size_t nPixels = dim[0] * dim[1];

  const size_t nInputObjectPDArrays = inputObjectPD->GetNumberOfArrays();
  const size_t nInputObjectCDArrays = inputObjectCD->GetNumberOfArrays();

  int status = 0;

  // Map Point Data
  for(size_t j = 0; j < nInputObjectPDArrays; j++) {
    auto inputArray = inputObjectPD->GetArray(j);
    auto outputArray
      = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
    outputArray->SetName(inputArray->GetName());
    outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
    outputArray->SetNumberOfTuples(nPixels);

    outputImagePD->AddArray(outputArray);

    switch(outputArray->GetDataType()) {
      vtkTemplateMacro(status = renderer->interpolateArray(
                         (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),

                         primitiveIdArray, barycentricCoordinates,
                         inputObjectConnectivityList,

                         (const VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
                         nPixels, inputArray->GetNumberOfComponents()));
    }

    if(!status)
      return 0;
  }

  // Map Cell Data
  for(size_t j = 0; j < nInputObjectCDArrays; j++) {
    auto inputArray = inputObjectCD->GetArray(j);
    auto outputArray
      = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
    outputArray->SetName(inputArray->GetName());
    outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
    outputArray->SetNumberOfTuples(nPixels);

    outputImagePD->AddArray(outputArray);

    switch(outputArray->GetDataType()) {
      vtkTemplateMacro(status = renderer->lookupArray(
                         (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),

                         primitiveIdArray,
                         (const VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
                         nPixels, inputArray->GetNumberOfComponents()));
    }

    if(!status)
      return 0;
  }

  return 1;
};
