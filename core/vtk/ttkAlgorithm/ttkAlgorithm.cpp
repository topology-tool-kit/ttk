#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <OrderDisambiguation.h>
#include <Triangulation.h>
#include <ttkTriangulationFactory.h>

#include <vtkCellTypes.h>
#include <vtkCommand.h>
#include <vtkDataSet.h>
#if TTK_ENABLE_MPI
#include <vtkCellData.h>
#include <vtkGenerateGlobalIds.h>
#include <vtkGhostCellsGenerator.h>
#endif
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationIntegerKey.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCompositeDataPipeline.h>

// Pass input type information key
#include <vtkInformationKey.h>
vtkInformationKeyMacro(ttkAlgorithm, SAME_DATA_TYPE_AS_INPUT_PORT, Integer);

// Constructor / Destructor
vtkStandardNewMacro(ttkAlgorithm);
ttkAlgorithm::ttkAlgorithm() = default;
ttkAlgorithm::~ttkAlgorithm() = default;

ttk::Triangulation *ttkAlgorithm::GetTriangulation(vtkDataSet *dataSet) {

  this->printMsg("Requesting triangulation for '"
                   + std::string(dataSet->GetClassName()) + "'",
                 ttk::debug::Priority::DETAIL);
#if TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    this->MPIPipelinePreconditioning(dataSet);
  }
#endif
  auto triangulation = ttkTriangulationFactory::GetTriangulation(
    this->debugLevel_, this->CompactTriangulationCacheSize, dataSet);
#if TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    if(triangulation) {
      this->MPITriangulationPreconditioning(triangulation, dataSet);
    }
  }
#endif
  if(triangulation)
    return triangulation;

  this->printErr("Unable to retrieve/initialize triangulation for '"
                 + std::string(dataSet->GetClassName()) + "'");

  return nullptr;
}

vtkDataArray *ttkAlgorithm::GetOptionalArray(const bool &enforceArrayIndex,
                                             const int &arrayIndex,
                                             const std::string &arrayName,
                                             vtkDataSet *const inputData,
                                             const int &inputPort) {

  vtkDataArray *optionalArray = nullptr;

  if(enforceArrayIndex)
    optionalArray = this->GetInputArrayToProcess(arrayIndex, inputData);

  if(!optionalArray) {
    this->SetInputArrayToProcess(arrayIndex, inputPort, 0, 0, arrayName.data());
    optionalArray = this->GetInputArrayToProcess(arrayIndex, inputData);
  }
  return optionalArray;
}

std::string ttkAlgorithm::GetOrderArrayName(vtkDataArray *const array) {
  return std::string(array->GetName()) + "_Order";
}

vtkDataArray *ttkAlgorithm::GetOrderArray(vtkDataSet *const inputData,
                                          const int scalarArrayIdx,
                                          const int orderArrayIdx,
                                          const bool enforceOrderArrayIdx) {

  auto isValidOrderArray = [](vtkDataArray *const array) {
    if(!array)
      return -4;

    if(array->GetNumberOfComponents() != 1)
      return -3;

    auto temp = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    if(array->GetDataType() != temp->GetDataType())
      return -2;

    const std::string name(array->GetName());
    if(name.size() < 6 || (name.rfind("_Order") != (name.size() - 6)))
      return -1;

    return 1;
  };

  if(enforceOrderArrayIdx) {
    auto orderArray = this->GetInputArrayToProcess(orderArrayIdx, inputData);
    switch(isValidOrderArray(orderArray)) {
      case -4: {
        this->printErr("Unable to retrieve enforced order array at idx "
                       + std::to_string(orderArrayIdx) + ".");
        return nullptr;
      }
      case -3: {
        this->printErr("Retrieved enforced order array `"
                       + std::string(orderArray->GetName())
                       + "` has more than one component.");
        return nullptr;
      }
      case -2: {
        this->printErr("Enforced order array `"
                       + std::string(orderArray->GetName())
                       + "` is of incorrect type.");
        auto temp = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
        this->printErr(" -> use `ttkArrayEditor` to convert data type to `"
                       + std::string(temp->GetDataTypeAsString()) + "`.");
        return nullptr;
      }
      default: {
        this->printMsg("Retrieved enforced order array `"
                         + std::string(orderArray->GetName()) + "`.",
                       ttk::debug::Priority::DETAIL);
        return orderArray;
      }
    }
  }

  auto scalarArray = this->GetInputArrayToProcess(scalarArrayIdx, inputData);
  if(!scalarArray) {
    this->printErr("Unable to retrieve input scalar array for idx "
                   + std::to_string(scalarArrayIdx) + ".");
    return nullptr;
  } else if(isValidOrderArray(scalarArray) == 1) {
    this->printMsg("Retrieved scalar array `"
                     + std::string(scalarArray->GetName())
                     + "` is already an order array.",
                   ttk::debug::Priority::DETAIL);
    return scalarArray;
  }

  auto orderArray = inputData
                      ->GetAttributesAsFieldData(this->GetInputArrayAssociation(
                        scalarArrayIdx, inputData))
                      ->GetArray(this->GetOrderArrayName(scalarArray).data());

  switch(isValidOrderArray(orderArray)) {
    case -4: {
      ttk::Timer timer;
      this->printWrn("No pre-existing order for array:");
      this->printWrn("  `" + std::string(scalarArray->GetName()) + "`.");

      this->printMsg("Initializing order array.", 0, 0, this->threadNumber_,
                     ttk::debug::LineMode::REPLACE);

      auto nVertices = scalarArray->GetNumberOfTuples();
      auto newOrderArray = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
      newOrderArray->SetName(this->GetOrderArrayName(scalarArray).data());
      newOrderArray->SetNumberOfComponents(1);
      newOrderArray->SetNumberOfTuples(nVertices);
#if TTK_ENABLE_MPI
      if(ttk::isRunningWithMPI()) {
        this->MPIPipelinePreconditioning(inputData);
        auto vtkGlobalPointIds = inputData->GetPointData()->GetGlobalIds();
        auto rankArray = inputData->GetPointData()->GetArray("RankArray");
        ttkTypeMacroAI(
          scalarArray->GetDataType(), vtkGlobalPointIds->GetDataType(),
          (ttk::produceOrdering<T0, T1>(
            ttkUtils::GetPointer<ttk::SimplexId>(newOrderArray),
            ttkUtils::GetPointer<T0>(scalarArray),
            ttkUtils::GetPointer<T1>(vtkGlobalPointIds),
            ttkUtils::GetPointer<int>(rankArray), nVertices, 500)));
      } else {
        switch(scalarArray->GetDataType()) {
          vtkTemplateMacro(ttk::preconditionOrderArray(
            nVertices,
            static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)),
            static_cast<ttk::SimplexId *>(
              ttkUtils::GetVoidPointer(newOrderArray)),
            this->threadNumber_));
        }
      }

#else
      switch(scalarArray->GetDataType()) {
        vtkTemplateMacro(ttk::preconditionOrderArray(
          nVertices,
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)),
          static_cast<ttk::SimplexId *>(
            ttkUtils::GetVoidPointer(newOrderArray)),
          this->threadNumber_));
      }
#endif

      // append order array temporarily to input
      inputData
        ->GetAttributesAsFieldData(
          this->GetInputArrayAssociation(scalarArrayIdx, inputData))
        ->AddArray(newOrderArray);
      this->printMsg("Initializing order array.", 1, timer.getElapsedTime(),
                     this->threadNumber_);

      this->printWrn("TIP: run `ttkArrayPreconditioning` first");
      this->printWrn("for improved performances :)");

      return newOrderArray;
    }

    case -3: {
      this->printErr(
        "Retrieved order array `" + std::string(orderArray->GetName())
        + "` for scalar array `" + std::string(scalarArray->GetName())
        + "` has more than one component.");
      return nullptr;
    }

    case -2: {
      this->printErr(
        "Retrieved order array `" + std::string(orderArray->GetName())
        + "` for scalar array `" + std::string(scalarArray->GetName())
        + "` is of incorrect type.");
      auto temp = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
      this->printErr(" -> use `ttkArrayEditor` to convert data type to `"
                     + std::string(temp->GetDataTypeAsString()) + "`.");
      return nullptr;
    }

    default: {
      this->printMsg(
        "Retrieved order array `" + std::string(orderArray->GetName())
          + "` for scalar array `" + std::string(scalarArray->GetName()) + "`.",
        ttk::debug::Priority::DETAIL);
      return orderArray;
    }
  }
}

ttk::SimplexId *
  ttkAlgorithm::GetIdentifierArrayPtr(const bool &enforceArrayIndex,
                                      const int &arrayIndex,
                                      const std::string &arrayName,
                                      vtkDataSet *const inputData,
                                      std::vector<ttk::SimplexId> &spareStorage,
                                      const int inputPort,
                                      const bool printErr) {

  // fetch data array
  const auto array = this->GetOptionalArray(
    enforceArrayIndex, arrayIndex, arrayName, inputData, inputPort);
  if(array == nullptr) {
    if(printErr) {
      this->printErr("Could not find the requested identifiers array");
    }
    return {};
  }
  if(array->GetNumberOfComponents() != 1) {
    if(printErr) {
      this->printErr("Identifiers field must have only one component!");
    }
    return {};
  }

#ifndef TTK_ENABLE_64BIT_IDS
  if(array->GetDataType() == VTK_ID_TYPE
     || array->GetDataType() == VTK_LONG_LONG) {
    this->printMsg(
      "Converting identifiers field from vtkIdType to SimplexId...");
    const auto nItems = array->GetNumberOfTuples();

    // fills the vector with the content of the data array converted to
    // ttk::SimplexId
    spareStorage.resize(nItems);
    for(vtkIdType i = 0; i < nItems; ++i) {
      spareStorage[i] = static_cast<ttk::SimplexId>(array->GetTuple1(i));
    }

    // return a pointer to the vector internal buffer
    return spareStorage.data();
  }
#else
  TTK_FORCE_USE(spareStorage);
#endif

  // return a pointer to the data array internal buffer
  return static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(array));
}

template <class vtkDataType>
int prepOutput(vtkInformation *info, const std::string &className) {
  auto output = vtkDataObject::GetData(info);
  if(!output || !output->IsA(className.data())) {
    auto newOutput = vtkSmartPointer<vtkDataType>::New();
    info->Set(vtkDataObject::DATA_OBJECT(), newOutput);
  }
  return 1;
}

vtkDataSet *ttkAlgorithm::GetOutput() {
  return this->GetOutput(0);
}

vtkDataSet *ttkAlgorithm::GetOutput(int port) {
  return vtkDataSet::SafeDownCast(this->GetOutputDataObject(port));
}

void ttkAlgorithm::SetInputData(vtkDataSet *input) {
  this->SetInputData(0, input);
}

void ttkAlgorithm::SetInputData(int index, vtkDataSet *input) {
  this->SetInputDataInternal(index, input);
}

void ttkAlgorithm::AddInputData(vtkDataSet *input) {
  this->AddInputData(0, input);
}

void ttkAlgorithm::AddInputData(int index, vtkDataSet *input) {
  this->AddInputDataInternal(index, input);
}

int ttkAlgorithm::RequestDataObject(vtkInformation *ttkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // for each output
  for(int i = 0; i < this->GetNumberOfOutputPorts(); ++i) {
    auto outInfo = outputVector->GetInformationObject(i);
    if(!outInfo) {
      this->printErr("Unable to retrieve output vtkDataObject at port "
                     + std::to_string(i));
      return 0;
    }

    auto outputPortInfo = this->GetOutputPortInformation(i);

    // always request output type again for dynamic filter outputs
    if(!this->FillOutputPortInformation(i, outputPortInfo)) {
      this->printErr("Unable to fill output port information at port "
                     + std::to_string(i));
      return 0;
    }

    if(outputPortInfo->Has(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT())) {
      // Set output data type to input data type at specified port
      auto inPortIndex
        = outputPortInfo->Get(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT());
      if(inPortIndex < 0 || inPortIndex >= this->GetNumberOfInputPorts()) {
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

#if TTK_ENABLE_MPI
void ttkAlgorithm::MPIPipelinePreconditioning(vtkDataSet *input) {

  vtkNew<vtkGenerateGlobalIds> globalIds;
  // If the global point id array doesn't exist, it is created
  if(input->GetPointData()->GetGlobalIds() == nullptr) {
    printWrn("Global ids haven't been produced in sequential, the parallel "
             "result may be different");
    globalIds->SetInputData(input);
    globalIds->Update();
    input->ShallowCopy(globalIds->GetOutputDataObject(0));
  }

  vtkNew<vtkGhostCellsGenerator> generator;
  if(!input->HasAnyGhostCells()) {
    generator->SetInputData(input);
    generator->BuildIfRequiredOff();
    generator->SetNumberOfGhostLayers(2);
    generator->Update();
    input->ShallowCopy(generator->GetOutputDataObject(0));
    input->GetPointData()->AddArray(
      generator->GetOutputDataObject(0)->GetGhostArray(0));
    input->GetCellData()->AddArray(
      generator->GetOutputDataObject(0)->GetGhostArray(1));
  }

  // If the RankArray array doesn't exist for pointdata, it is created
  if(input->GetPointData()->GetArray("RankArray") == nullptr) {
    int vertexNumber = input->GetNumberOfPoints();
    std::vector<int> rankArray(vertexNumber, 0);
    double *boundingBox = input->GetBounds();
    ttk::produceRankArray(rankArray,
                          static_cast<long int *>(ttkUtils::GetVoidPointer(
                            input->GetPointData()->GetGlobalIds())),
                          static_cast<unsigned char *>(ttkUtils::GetVoidPointer(
                            input->GetPointData()->GetArray("vtkGhostType"))),
                          vertexNumber, boundingBox);

    vtkNew<vtkIntArray> vtkRankArray{};
    vtkRankArray->SetName("RankArray");
    vtkRankArray->SetNumberOfComponents(1);
    vtkRankArray->SetNumberOfTuples(vertexNumber);

    for(int i = 0; i < vertexNumber; i++) {
      vtkRankArray->SetComponent(i, 0, rankArray[i]);
    }

    input->GetPointData()->AddArray(vtkRankArray);
  }

  // If the RankArray array doesn't exist for celldata, it is created
  if(input->GetCellData()->GetArray("RankArray") == nullptr) {
    int cellNumber = input->GetNumberOfCells();
    std::vector<int> cellsRankArray(cellNumber, 0);
    double *boundingBox = input->GetBounds();
    ttk::produceRankArray(cellsRankArray,
                          static_cast<long int *>(ttkUtils::GetVoidPointer(
                            input->GetCellData()->GetGlobalIds())),
                          static_cast<unsigned char *>(ttkUtils::GetVoidPointer(
                            input->GetCellData()->GetArray("vtkGhostType"))),
                          cellNumber, boundingBox);

    vtkNew<vtkIntArray> vtkCellsRankArray{};
    vtkCellsRankArray->SetName("RankArray");
    vtkCellsRankArray->SetNumberOfComponents(1);
    vtkCellsRankArray->SetNumberOfTuples(cellNumber);

    for(int i = 0; i < cellNumber; i++) {
      vtkCellsRankArray->SetComponent(i, 0, cellsRankArray[i]);
    }

    input->GetCellData()->AddArray(vtkCellsRankArray);
  }
}

void ttkAlgorithm::MPITriangulationPreconditioning(
  ttk::Triangulation *triangulation, vtkDataSet *input) {
  triangulation->setGlobalIdsArray(static_cast<long int *>(
    ttkUtils::GetVoidPointer(input->GetPointData()->GetGlobalIds())));
  triangulation->setGlobalIdsArray(static_cast<long int *>(
    ttkUtils::GetVoidPointer(input->GetPointData()->GetGlobalIds())));
  triangulation->preconditionDistributedVertices();
  triangulation->setRankArray(static_cast<int *>(
    ttkUtils::GetVoidPointer(input->GetPointData()->GetArray("RankArray"))));
  // provide "GlobalCellIds" & "vtkGhostType" cell data array to
  // the triangulation
  const auto cd{input->GetCellData()};
  if(cd == nullptr) {
    triangulation->printWrn("No cell data on input object");
  }
  if(cd != nullptr) {
    triangulation->setGlobalIds(
      ttkUtils::GetPointer<ttk::LongSimplexId>(cd->GetGlobalIds()),
      ttkUtils::GetPointer<unsigned char>(
        cd->GetArray(vtkCellData::GhostArrayName())));
  }
}
#endif
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
