#include <ttkConnectedComponents.h>

#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vtkFloatArray.h>
#include <vtkIntArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkConnectedComponents);

ttkConnectedComponents::ttkConnectedComponents() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);

  // Suppress warning if one does not set the optional input array
  this->SetInputArrayToProcess(0, 0, 0, 0, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
}

int ttkConnectedComponents::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkConnectedComponents::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkConnectedComponents::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  // Fetch Input Data
  auto inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  const auto nPoints = inputDataSet->GetNumberOfPoints();

  auto featureMask = this->GetInputArrayToProcess(0, inputVector);

  if(featureMask && this->GetInputArrayAssociation(0, inputVector) != 0)
    return !this->printErr("Input array needs to be a point data array.");

  if(featureMask && featureMask->GetNumberOfComponents() != 1)
    return !this->printErr("Input array needs to be a scalar array.");

  auto triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  // Allocate Output Label Data
  auto componentIds = vtkSmartPointer<vtkIntArray>::New();
  componentIds->SetName("ComponentId");
  componentIds->SetNumberOfTuples(nPoints);

  auto componentSize = vtkSmartPointer<vtkFloatArray>::New();
  componentSize->SetName("ComponentSize");
  if(this->AugmentSegmentationWithComponentSize)
    componentSize->SetNumberOfTuples(nPoints);

  // Compute Connected Components
  std::vector<ttk::ConnectedComponents::Component> components;
  {
    int status = 0;

    // initialize component ids
    ttkTypeMacroA(featureMask ? featureMask->GetDataType() : VTK_INT,
                  (status = this->initializeComponentIds<T0>(
                     ttkUtils::GetPointer<int>(componentIds), nPoints,
                     ttkUtils::GetPointer<const T0>(featureMask),
                     this->BackgroundThreshold)));
    if(!status)
      return 0;

    // compute connected components (componentIds now store component idx)
    this->preconditionTriangulation(triangulation);
    ttkTypeMacroT(triangulation->getType(),
                  (status = this->computeConnectedComponents<T0>(
                     components, ttkUtils::GetPointer<int>(componentIds),
                     static_cast<const T0 *>(triangulation->getData()))));
    if(!status)
      return 0;
  }

  // Augment Segmentation
  {
    int status = 0;

    // map component size to segmentation
    if(this->AugmentSegmentationWithComponentSize) {
      status = this->mapData(
        components, ttkUtils::GetPointer<int>(componentIds), nPoints,
        [](const std::vector<ttk::ConnectedComponents::Component> &components_,
           const int cIdx) { return cIdx < 0 ? 0.0 : components_[cIdx].size; },
        ttkUtils::GetPointer<float>(componentSize));
      if(!status)
        return 0;
    }

    // map component idx to component id
    status = this->mapData(
      components, ttkUtils::GetPointer<int>(componentIds), nPoints,
      [](const std::vector<ttk::ConnectedComponents::Component> &components_,
         const int cIdx) { return cIdx < 0 ? -1 : components_[cIdx].id; },
      ttkUtils::GetPointer<int>(componentIds));
    if(!status)
      return 0;
  }

  // Segmentation Output
  {
    auto outputDataSet = vtkDataSet::GetData(outputVector, 0);
    outputDataSet->ShallowCopy(inputDataSet);
    outputDataSet->GetPointData()->AddArray(componentIds);
    if(this->AugmentSegmentationWithComponentSize)
      outputDataSet->GetPointData()->AddArray(componentSize);
  }

  // Components Output
  {
    const int nComponents = components.size();
    auto outputComponents = vtkPolyData::GetData(outputVector, 1);

    // points
    {
      auto sizeArray = vtkSmartPointer<vtkFloatArray>::New();
      sizeArray->SetName("Size");
      sizeArray->SetNumberOfTuples(nComponents);
      auto sizeArrayData = ttkUtils::GetPointer<float>(sizeArray);

      auto idArray = vtkSmartPointer<vtkIntArray>::New();
      idArray->SetName("ComponentId");
      idArray->SetNumberOfTuples(nComponents);
      auto idArrayData = ttkUtils::GetPointer<int>(idArray);

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetDataTypeToFloat();
      points->SetNumberOfPoints(nComponents);
      auto pointsData = ttkUtils::GetPointer<float>(points->GetData());
      for(int i = 0, j = 0; i < nComponents; i++) {
        const auto &c = components[i];
        pointsData[j++] = c.center[0];
        pointsData[j++] = c.center[1];
        pointsData[j++] = c.center[2];

        sizeArrayData[i] = c.size;
        idArrayData[i] = c.id;
      }

      outputComponents->SetPoints(points);
      auto pd = outputComponents->GetPointData();
      pd->AddArray(sizeArray);
      pd->AddArray(idArray);
    }

    // cells
    {
      auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
      connectivityArray->SetNumberOfTuples(nComponents);
      auto connectivityArrayData = ttkUtils::GetPointer<int>(connectivityArray);
      for(int i = 0; i < nComponents; i++)
        connectivityArrayData[i] = i;

      auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
      offsetArray->SetNumberOfTuples(nComponents + 1);
      auto offsetArrayData = ttkUtils::GetPointer<int>(offsetArray);
      for(int i = 0; i <= nComponents; i++)
        offsetArrayData[i] = i;

      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetData(offsetArray, connectivityArray);

      outputComponents->SetVerts(cellArray);
    }

    // Copy Field Data
    outputComponents->GetFieldData()->ShallowCopy(inputDataSet->GetFieldData());
  }

  // return success
  return 1;
}
