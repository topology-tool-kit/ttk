#include <ttkAlgorithm.h>
#include <ttkUtils.h>

#include <vtkDataObject.h> // For output port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkSmartPointer.h>

#include <vtkInformationVector.h>

#include <vtkCompositeDataPipeline.h>

#include <vtkFieldData.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// TODO: use a class here to add semantic about the four fields
// and clear access methods
typedef std::unordered_map<void *,
                           std::tuple<ttk::Triangulation,
                                      vtkObject *,
                                      vtkSmartPointer<vtkCommand>,
                                      vtkMTimeType>>
  DataSetToTriangulationMapType;
DataSetToTriangulationMapType ttkAlgorithm::DataSetToTriangulationMap;

struct ttkOnDeleteCommand : public vtkCommand {
  bool deleteEventFired{false};
  vtkObject *owner_;
  DataSetToTriangulationMapType *dataSetToTriangulationMap_;

  static ttkOnDeleteCommand *New() {
    return new ttkOnDeleteCommand;
  }
  vtkTypeMacro(ttkOnDeleteCommand, vtkCommand);

  void Init(vtkObject *owner,
            DataSetToTriangulationMapType *dataSetToTriangulationMap) {
    this->owner_ = owner;
    this->owner_->AddObserver(vtkCommand::DeleteEvent, this, 1);
    this->dataSetToTriangulationMap_ = dataSetToTriangulationMap;
  }
  ~ttkOnDeleteCommand() {
    if(!this->deleteEventFired)
      this->owner_->RemoveObserver(this);
  }

  void Execute(vtkObject *, unsigned long eventId, void *callData) override {
    this->deleteEventFired = true;

    void *key = this->owner_;
    if(this->owner_->IsA("vtkImageData"))
      key = vtkImageData::SafeDownCast(owner_)->GetScalarPointer();

    auto it = this->dataSetToTriangulationMap_->find(key);
    if(it != this->dataSetToTriangulationMap_->end())
      this->dataSetToTriangulationMap_->erase(it);
  }
};

// Pass input type information key
#include <vtkInformationKey.h>
vtkInformationKeyMacro(ttkAlgorithm, SAME_DATA_TYPE_AS_INPUT_PORT, Integer);

vtkStandardNewMacro(ttkAlgorithm);

ttkAlgorithm::ttkAlgorithm() {
}
ttkAlgorithm::~ttkAlgorithm() {
}

ttk::Triangulation *ttkAlgorithm::FindTriangulation(void *key) {
  auto it = ttkAlgorithm::DataSetToTriangulationMap.find(key);
  if(it != ttkAlgorithm::DataSetToTriangulationMap.end()) {
    auto triangulation = &std::get<0>(it->second);
    auto owner = std::get<1>(it->second);
    bool valid = false;

    // check if triangulation is still valid
    if(owner->IsA("vtkImageData")) {
      auto ownerAsID = vtkImageData::SafeDownCast(owner);

      std::vector<int> ttkDimensions;
      triangulation->getGridDimensions(ttkDimensions);

      int vtkDimensions[3];
      ownerAsID->GetDimensions(vtkDimensions);

      valid = (vtkDimensions[0] == ttkDimensions[0])
              && (vtkDimensions[1] == ttkDimensions[1])
              && (vtkDimensions[2] == ttkDimensions[2]);
    } else if(owner->IsA("vtkCellArray")) {
      valid = std::get<3>(it->second) == owner->GetMTime();
    }

    if(valid) {
      this->printMsg("Returning already initilized triangulation",
                     ttk::debug::Priority::DETAIL);
      triangulation->setDebugLevel(this->debugLevel_);
      return triangulation;
    } else {
      this->printMsg(
        "Chached triangulation no longer valid", ttk::debug::Priority::DETAIL);
      ttkAlgorithm::DataSetToTriangulationMap.erase(it);
    }
  }

  return nullptr;
}

ttk::Triangulation *ttkAlgorithm::InitTriangulation(void *key,
                                                    vtkObject *owner,
                                                    vtkPoints *points,
                                                    vtkCellArray *cells) {
  ttkAlgorithm::DataSetToTriangulationMap.insert(
    {key,
     {ttk::Triangulation(), owner, vtkSmartPointer<ttkOnDeleteCommand>::New(),
      owner->GetMTime()}});

  auto it = ttkAlgorithm::DataSetToTriangulationMap.find(key);
  auto triangulation = &std::get<0>(it->second);
  triangulation->setDebugLevel(this->debugLevel_);

  // Delete callback
  {
    ttkOnDeleteCommand::SafeDownCast(std::get<2>(it->second))
      ->Init(owner, &ttkAlgorithm::DataSetToTriangulationMap);
  }

  // Initialize Explicit Triangulation
  if(points && cells) {
    this->printMsg("Initializing Explicit Triangulation", 0,
                   ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);

    // Points
    {
      auto pointDataType = points->GetDataType();
      if(pointDataType != VTK_FLOAT && pointDataType != VTK_DOUBLE) {
        this->printErr("Unable to initialize 'ttk::Triangulation' for point "
                       "precision other than 'float' or 'double'.");
        ttkAlgorithm::DataSetToTriangulationMap.erase(it);
        return nullptr;
      }

      void *pointDataArray = ttkUtils::GetVoidPointer(points);
      triangulation->setInputPoints(points->GetNumberOfPoints(), pointDataArray,
                                    pointDataType == VTK_DOUBLE);
    }

    // Cells
    if(cells->GetNumberOfCells() > 0) {
#ifdef TTK_CELL_ARRAY_NEW
      auto connectivity
        = (vtkIdType *)ttkUtils::GetVoidPointer(cells->GetConnectivityArray());
      auto offsets
        = (vtkIdType *)ttkUtils::GetVoidPointer(cells->GetOffsetsArray());
      triangulation->setInputCells(
        cells->GetNumberOfCells(), connectivity, offsets);
#else
      triangulation->setInputCells(
        cells->GetNumberOfCells(), cells->GetData()->GetPointer(0));
#endif
    }

    this->printMsg(
      "Initializing Explicit Triangulation", 1, ttk::debug::Priority::DETAIL);
  } else {
    // Initialize Implicit Triangulation
    this->printMsg("Initializing Implicit Triangulation", 0,
                   ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);

    auto ownerAsID = vtkImageData::SafeDownCast(owner);

    int extents[6];
    ownerAsID->GetExtent(extents);

    double origin[3];
    ownerAsID->GetOrigin(origin);

    double spacing[3];
    ownerAsID->GetSpacing(spacing);

    int gridDimensions[3];
    ownerAsID->GetDimensions(gridDimensions);

    double firstPoint[3];
    firstPoint[0] = origin[0] + extents[0] * spacing[0];
    firstPoint[1] = origin[1] + extents[2] * spacing[1];
    firstPoint[2] = origin[2] + extents[4] * spacing[2];

    triangulation->setInputGrid(
      firstPoint[0], firstPoint[1], firstPoint[2], spacing[0], spacing[1],
      spacing[2], gridDimensions[0], gridDimensions[1], gridDimensions[2]);

    this->printMsg(
      "Initializing Implicit Triangulation", 1, ttk::debug::Priority::DETAIL);
  }

  this->printMsg(
    "Number of registered triangulations: "
      + std::to_string(ttkAlgorithm::DataSetToTriangulationMap.size()),
    ttk::debug::Priority::VERBOSE);

  return triangulation;
}

ttk::Triangulation *ttkAlgorithm::GetTriangulation(vtkDataSet *dataSet) {

  this->printMsg("Requesting triangulation for '"
                   + std::string(dataSet->GetClassName()) + "'",
                 ttk::debug::Priority::DETAIL);

  switch(dataSet->GetDataObjectType()) {
    // =========================================================================
    case VTK_UNSTRUCTURED_GRID: {
      auto dataSetAsUG = vtkUnstructuredGrid::SafeDownCast(dataSet);
      auto points = dataSetAsUG->GetPoints();
      auto cells = dataSetAsUG->GetCells();

      // check if triangulation already exists or has to be updated
      ttk::Triangulation *triangulation = this->FindTriangulation(cells);

      // otherwise create new triangulation
      if(!triangulation)
        triangulation = this->InitTriangulation(cells, cells, points, cells);

      // return triangulation
      return triangulation;
    }

    // =========================================================================
    case VTK_POLY_DATA: {
      auto dataSetAsPD = vtkPolyData::SafeDownCast(dataSet);
      auto points = dataSetAsPD->GetPoints();
      auto polyCells = dataSetAsPD->GetPolys();
      auto lineCells = dataSetAsPD->GetLines();

      // check if triangulation for polycells already exists or has to be
      // updated
      ttk::Triangulation *triangulation = this->FindTriangulation(polyCells);

      // check if triangulation for lineCells already exists or has to be
      // updated
      if(!triangulation)
        triangulation = this->FindTriangulation(lineCells);

      // if not create triangulation
      if(!triangulation) {
        vtkCellArray *cellArray = nullptr;
        if(polyCells->GetNumberOfCells() > 0)
          cellArray = polyCells;
        else if(lineCells->GetNumberOfCells() > 0)
          cellArray = lineCells;

        if(cellArray != nullptr) {
          triangulation
            = this->InitTriangulation(cellArray, cellArray, points, cellArray);
        }
      }

      if(!triangulation) {
        this->printErr("Unable to initialize triangulation for vtkPolyData "
                       "without any cells.");
      }

      return triangulation;
    }

    // =========================================================================
    case VTK_IMAGE_DATA: {
      auto dataSetAsID = vtkImageData::SafeDownCast(dataSet);

      // check if triangulation already exists or has to be updated
      ttk::Triangulation *triangulation
        = this->FindTriangulation(dataSetAsID->GetScalarPointer());

      // otherwise initialize triangulation
      if(!triangulation) {
        triangulation
          = this->InitTriangulation(dataSetAsID->GetScalarPointer(), dataSet);
      }

      // return triangulation
      return triangulation;
    }

    // UNSUPPORTED DATA TYPE
    // =========================================================================
    default: {
    }
  }

  this->printErr("Unable to retrieve/initialize triangulation for '"
                 + std::string(dataSet->GetClassName()) + "'");

  return nullptr;
}

template <class vtkDataType>
int prepOutput(vtkInformation *info, std::string className) {
  auto output = vtkDataObject::GetData(info);
  if(!output || !output->IsA(className.data())) {
    auto newOutput = vtkSmartPointer<vtkDataType>::New();
    info->Set(vtkDataObject::DATA_OBJECT(), newOutput);
  }
  return 1;
}

int ttkAlgorithm::RequestDataObject(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // for each output
  for(int i = 0; i < this->GetNumberOfOutputPorts(); ++i) {
    auto outInfo = outputVector->GetInformationObject(i);

    auto outputPortInfo = this->GetOutputPortInformation(i);

    // always request output type again for dynamic filter outputs
    this->FillOutputPortInformation(i, outputPortInfo);

    if(outputPortInfo->Has(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT())) {
      // Set output data type to input data type at specified port
      auto inPortIndex
        = outputPortInfo->Get(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT());
      if(inPortIndex >= this->GetNumberOfInputPorts()) {
        this->printErr("Input port index " + std::to_string(inPortIndex)
                       + " specified by 'SAME_DATA_TYPE_AS_INPUT_PORT' key of "
                         "output port is out of range ("
                       + std::to_string(this->GetNumberOfInputPorts())
                       + " input ports).");
        return 0;
      }
      auto inInfo = inputVector[inPortIndex]->GetInformationObject(0);
      if(!inInfo) {
        this->printErr(
          "No information object at port " + std::to_string(inPortIndex)
          + " specified by 'SAME_DATA_TYPE_AS_INPUT_PORT' key of output port.");
        return 0;
      }

      auto input = vtkDataObject::GetData(inInfo);
      auto output = vtkDataObject::GetData(outInfo);

      if(!output || !output->IsA(input->GetClassName())) {
        auto newOutput
          = vtkSmartPointer<vtkDataObject>::Take(input->NewInstance());
        outputPortInfo->Set(
          vtkDataObject::DATA_TYPE_NAME(), input->GetClassName());
        outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutput);
      }
    } else {
      // Explicitly create output by data type name
      if(!outputPortInfo->Has(vtkDataObject::DATA_TYPE_NAME())) {
        this->printErr("DATA_TYPE_NAME of output port " + std::to_string(i)
                       + " not specified");
        return 0;
      }
      std::string outputType
        = outputPortInfo->Get(vtkDataObject::DATA_TYPE_NAME());

      if(outputType == "vtkUnstructuredGrid") {
        prepOutput<vtkUnstructuredGrid>(outInfo, outputType);
      } else if(outputType == "vtkPolyData") {
        prepOutput<vtkPolyData>(outInfo, outputType);
      } else if(outputType == "vtkMultiBlockDataSet") {
        prepOutput<vtkMultiBlockDataSet>(outInfo, outputType);
      } else if(outputType == "vtkTable") {
        prepOutput<vtkTable>(outInfo, outputType);
      } else if(outputType == "vtkImageData") {
        prepOutput<vtkImageData>(outInfo, outputType);
      } else {
        this->printErr("Unsupported data type for output[" + std::to_string(i)
                       + "]: " + outputType);
        return 0;
      }
    }

    this->printMsg(
      "Created '"
        + std::string(outputPortInfo->Get(vtkDataObject::DATA_TYPE_NAME()))
        + "' at output port " + std::to_string(i),
      ttk::debug::Priority::VERBOSE);
  }

  return 1;
}

//==============================================================================
int ttkAlgorithm::ProcessRequest(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
  // 1. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA_OBJECT())) {
    this->printMsg(
      "Processing REQUEST_DATA_OBJECT", ttk::debug::Priority::VERBOSE);
    return this->RequestDataObject(request, inputVector, outputVector);
  }

  // 2. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_INFORMATION())) {
    this->printMsg(
      "Processing REQUEST_INFORMATION", ttk::debug::Priority::VERBOSE);
    return this->RequestInformation(request, inputVector, outputVector);
  }

  // 3. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_UPDATE_TIME())) {
    this->printMsg(
      "Processing REQUEST_UPDATE_TIME", ttk::debug::Priority::VERBOSE);
    return this->RequestUpdateTime(request, inputVector, outputVector);
  }

  // 4. Pass
  if(request->Has(
       vtkCompositeDataPipeline::REQUEST_TIME_DEPENDENT_INFORMATION())) {
    this->printMsg("Processing REQUEST_TIME_DEPENDENT_INFORMATION",
                   ttk::debug::Priority::VERBOSE);
    return this->RequestUpdateTimeDependentInformation(
      request, inputVector, outputVector);
  }

  // 5. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_UPDATE_EXTENT())) {
    this->printMsg(
      "Processing REQUEST_UPDATE_EXTENT", ttk::debug::Priority::VERBOSE);
    return this->RequestUpdateExtent(request, inputVector, outputVector);
  }

  // 6. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA_NOT_GENERATED())) {
    this->printMsg(
      "Processing REQUEST_DATA_NOT_GENERATED", ttk::debug::Priority::VERBOSE);
    return this->RequestDataNotGenerated(request, inputVector, outputVector);
  }

  // 7. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA())) {
    this->printMsg("Processing REQUEST_DATA", ttk::debug::Priority::VERBOSE);
    this->printMsg(ttk::debug::Separator::L0);
    return this->RequestData(request, inputVector, outputVector);
  }

  this->printErr("Unsupported pipeline pass:");
  request->Print(cout);

  return 0;
}
