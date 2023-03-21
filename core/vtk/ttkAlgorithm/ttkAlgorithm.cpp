#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <OrderDisambiguation.h>
#include <Triangulation.h>
#include <ttkTriangulationFactory.h>

#include <vtkCellTypes.h>
#include <vtkCommand.h>
#include <vtkDataSet.h>

#ifdef TTK_ENABLE_MPI
#include <Identifiers.h>
#include <vtkCellData.h>
#include <vtkGhostCellsGenerator.h>
#endif // TTK_ENABLE_MPI

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
#ifdef TTK_ENABLE_MPI
  if(ttk::hasInitializedMPI()) {
    if(!hasMPISupport_) {
      printErr(
        "MPI is not supported for this filter, the results will be incorrect");
    }
    this->MPIGhostPipelinePreconditioning(dataSet);
  }
#endif // TTK_ENABLE_MPI

  auto triangulation = ttkTriangulationFactory::GetTriangulation(
    this->debugLevel_, this->CompactTriangulationCacheSize, dataSet);

#ifdef TTK_ENABLE_MPI
  if(ttk::hasInitializedMPI()) {
    std::vector<int> tmp{};
    this->MPIPipelinePreconditioning(dataSet, tmp, triangulation);
    this->MPITriangulationPreconditioning(triangulation, dataSet);
  }
#endif // TTK_ENABLE_MPI

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
#ifdef TTK_ENABLE_MPI
      std::vector<int> neighbors;
      if(ttk::hasInitializedMPI()) {
        this->MPIGhostPipelinePreconditioning(inputData);
        this->MPIPipelinePreconditioning(inputData, neighbors, nullptr);
      }
      if(ttk::isRunningWithMPI()) {
        const auto triangulation{this->GetTriangulation(inputData)};
        ttkTypeMacroA(scalarArray->GetDataType(),
                      (ttk::produceOrdering<T0>(
                        ttkUtils::GetPointer<ttk::SimplexId>(newOrderArray),
                        ttkUtils::GetPointer<T0>(scalarArray),
                        [triangulation](const ttk::SimplexId a) {
                          return triangulation->getVertexGlobalId(a);
                        },
                        [triangulation](const ttk::SimplexId a) {
                          return triangulation->getVertexRank(a);
                        },
                        [triangulation](const ttk::SimplexId a) {
                          return triangulation->getVertexLocalId(a);
                        },
                        nVertices, 500, neighbors)));
      } else
#endif // TTK_ENABLE_MPI

      {
        switch(scalarArray->GetDataType()) {
          vtkTemplateMacro(ttk::preconditionOrderArray(
            nVertices,
            static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)),
            static_cast<ttk::SimplexId *>(
              ttkUtils::GetVoidPointer(newOrderArray)),
            this->threadNumber_));
        }
      }

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
#endif // TTK_ENABLE_64BIT_IDS

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

#ifdef TTK_ENABLE_MPI

int ttkAlgorithm::updateMPICommunicator(vtkDataSet *input) {
  int isEmpty
    = input->GetNumberOfCells() == 0 || input->GetNumberOfPoints() == 0;
  int oldSize = ttk::MPIsize_;
  int oldRank = ttk::MPIrank_;
  MPI_Comm_split(MPI_COMM_WORLD, isEmpty, 0, &ttk::MPIcomm_);
  MPI_Comm_rank(ttk::MPIcomm_, &ttk::MPIrank_);
  MPI_Comm_size(ttk::MPIcomm_, &ttk::MPIsize_);
  if(oldSize != ttk::MPIsize_) {
    std::vector<int> newToOldRanks(ttk::MPIsize_);
    MPI_Allgather(&oldRank, 1, MPI_INTEGER, newToOldRanks.data(), 1,
                  MPI_INTEGER, ttk::MPIcomm_);
    std::map<int, int> oldToNewRanks;
    for(int i = 0; i < ttk::MPIsize_; i++) {
      oldToNewRanks[newToOldRanks[i]] = i;
    }
    int *vertexRankArray
      = ttkUtils::GetPointer<int>(input->GetPointData()->GetArray("RankArray"));
    if(vertexRankArray != nullptr) {
      for(int i = 0; i < input->GetNumberOfPoints(); i++) {
        vertexRankArray[i] = oldToNewRanks[vertexRankArray[i]];
      }
    }
    int *cellRankArray
      = ttkUtils::GetPointer<int>(input->GetCellData()->GetArray("RankArray"));
    if(cellRankArray != nullptr) {
      for(int i = 0; i < input->GetNumberOfCells(); i++) {
        cellRankArray[i] = oldToNewRanks[cellRankArray[i]];
      }
    }
  }
  this->setDebugMsgPrefix(debugMsgNamePrefix_);
  return isEmpty;
}

bool ttkAlgorithm::checkGlobalIdValidity(ttk::LongSimplexId *globalIds,
                                         ttk::SimplexId simplexNumber,
                                         unsigned char *ghost,
                                         int *rankArray) {
  ttk::SimplexId ghostNumber = 0;
  if(rankArray != nullptr) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for reduction(+ : ghostNumber)
#endif // TTK_ENABLE_OPENMP
    for(ttk::SimplexId i = 0; i < simplexNumber; i++) {
      if(rankArray[i] != ttk::MPIrank_) {
        ghostNumber++;
      }
    }
  } else {
    if(ghost != nullptr) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for reduction(+ : ghostNumber)
#endif // TTK_ENABLE_OPENMP
      for(ttk::SimplexId i = 0; i < simplexNumber; i++) {
        if(ghost[i] == 1) {
          ghostNumber++;
        }
      }
    }
  }

  ttk::SimplexId realSimplexNumber = simplexNumber - ghostNumber;
  auto minmax = std::minmax_element(globalIds, globalIds + simplexNumber);
  ttk::LongSimplexId min = globalIds[minmax.first - globalIds];
  ttk::LongSimplexId max = globalIds[minmax.second - globalIds];
  ttk::SimplexId globalSimplexNumber;
  ttk::LongSimplexId globalMin;
  ttk::LongSimplexId globalMax;
  MPI_Allreduce(&realSimplexNumber, &globalSimplexNumber, 1,
                ttk::getMPIType(realSimplexNumber), MPI_SUM, ttk::MPIcomm_);
  MPI_Allreduce(
    &min, &globalMin, 1, ttk::getMPIType(min), MPI_MIN, ttk::MPIcomm_);
  MPI_Allreduce(
    &max, &globalMax, 1, ttk::getMPIType(max), MPI_MAX, ttk::MPIcomm_);

  return (globalSimplexNumber == globalMax + 1 && globalMin == 0);
};

int ttkAlgorithm::GenerateGlobalIds(
  vtkDataSet *input,
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> &vertGtoL,
  std::vector<int> &neighborRanks) {

  ttk::Identifiers identifiers;

  vtkNew<vtkIdTypeArray> vtkVertexIdentifiers{};

  vtkNew<vtkIdTypeArray> vtkCellIdentifiers{};
  vtkVertexIdentifiers->SetName("GlobalPointIds");
  vtkVertexIdentifiers->SetNumberOfComponents(1);
  vtkVertexIdentifiers->SetNumberOfTuples(input->GetNumberOfPoints());
  vtkVertexIdentifiers->Fill(-1);

  identifiers.setVertexIdentifiers(
    ttkUtils::GetPointer<ttk::LongSimplexId>(vtkVertexIdentifiers));

  vtkCellIdentifiers->SetName("GlobalCellIds");
  vtkCellIdentifiers->SetNumberOfComponents(1);
  vtkCellIdentifiers->SetNumberOfTuples(input->GetNumberOfCells());
  vtkCellIdentifiers->Fill(-1);

  identifiers.setCellIdentifiers(
    ttkUtils::GetPointer<ttk::LongSimplexId>(vtkCellIdentifiers));

  int vertexNumber = input->GetNumberOfPoints();
  identifiers.setVertexNumber(vertexNumber);
  int cellNumber = input->GetNumberOfCells();
  identifiers.setCellNumber(cellNumber);
  int status = 0;

  double *boundingBox = input->GetBounds();
  identifiers.setBounds(boundingBox);
  identifiers.initializeNeighbors(boundingBox, neighborRanks);
  if(ttk::isRunningWithMPI()) {
    switch(input->GetDataObjectType()) {
      case VTK_UNSTRUCTURED_GRID:
      case VTK_POLY_DATA: {

        identifiers.setOutdatedGlobalPointIds(
          ttkUtils::GetPointer<ttk::LongSimplexId>(
            input->GetPointData()->GetGlobalIds()));
        identifiers.setOutdatedGlobalCellIds(
          ttkUtils::GetPointer<ttk::LongSimplexId>(
            input->GetCellData()->GetGlobalIds()));
        identifiers.setVertexRankArray(ttkUtils::GetPointer<int>(
          input->GetPointData()->GetArray("RankArray")));
        int *cellRankArray = ttkUtils::GetPointer<int>(
          input->GetCellData()->GetArray("RankArray"));
        identifiers.setCellRankArray(cellRankArray);
        identifiers.setVertGhost(ttkUtils::GetPointer<unsigned char>(
          input->GetPointData()->GetArray("vtkGhostType")));
        unsigned char *cellGhost = ttkUtils::GetPointer<unsigned char>(
          input->GetCellData()->GetArray("vtkGhostType"));
        identifiers.setCellGhost(cellGhost);
        vtkPointSet *pointSet = vtkPointSet::SafeDownCast(input);
        identifiers.setPointSet(static_cast<float *>(
          ttkUtils::GetVoidPointer(pointSet->GetPoints())));
        vtkCellArray *cells = nullptr;
        switch(input->GetDataObjectType()) {
          case VTK_UNSTRUCTURED_GRID: {
            auto dataSetAsUG = vtkUnstructuredGrid::SafeDownCast(input);
            cells = dataSetAsUG->GetCells();
            break;
          }
          case VTK_POLY_DATA: {
            auto dataSetAsPD = vtkPolyData::SafeDownCast(input);
            cells
              = dataSetAsPD->GetNumberOfPolys() > 0   ? dataSetAsPD->GetPolys()
                : dataSetAsPD->GetNumberOfLines() > 0 ? dataSetAsPD->GetLines()
                                                      : dataSetAsPD->GetVerts();
            break;
          }
          default: {
            this->printErr("Unable to get cells for `"
                           + std::string(input->GetClassName()) + "`");
          }
        }
        if(cells == nullptr) {
          return 0;
        }
        if(!cells->IsStorage64Bit()) {
          if(cells->CanConvertTo64BitStorage()) {
            this->printWrn("Converting the cell array to 64-bit storage");
            bool success = cells->ConvertTo64BitStorage();
            if(!success) {
              this->printErr(
                "Error converting the provided cell array to 64-bit storage");
              return -1;
            }
          } else {
            this->printErr(
              "Cannot convert the provided cell array to 64-bit storage");
            return -1;
          }
        }

        identifiers.setConnectivity(ttkUtils::GetPointer<ttk::LongSimplexId>(
          cells->GetConnectivityArray()));

        std::vector<std::vector<ttk::SimplexId>> pointsToCells(vertexNumber);
        vtkIdList *cellList = vtkIdList::New();
        if(cellRankArray != nullptr) {
          for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
            input->GetPointCells(i, cellList);
            for(int j = 0; j < cellList->GetNumberOfIds(); j++) {
              if(cellRankArray[cellList->GetId(j)] == ttk::MPIrank_) {
                pointsToCells[i].push_back(cellList->GetId(j));
              }
            }
          }
        } else {
          for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
            input->GetPointCells(i, cellList);
            for(int j = 0; j < cellList->GetNumberOfIds(); j++) {
              if(cellGhost[cellList->GetId(j)] == 0) {
                pointsToCells[i].push_back(cellList->GetId(j));
              }
            }
          }
        }
        identifiers.setPointsToCells(pointsToCells);

        identifiers.initializeMPITypes();
        identifiers.setVertGtoL(&vertGtoL);
        vtkIdList *pointCell = vtkIdList::New();
        input->GetCellPoints(0, pointCell);
        int nbPoints = pointCell->GetNumberOfIds();
        identifiers.setDomainDimension(nbPoints - 1);
        identifiers.buildKDTree();
        status = identifiers.executePolyData();
        break;
      }
      case VTK_IMAGE_DATA: {
        vtkImageData *data = vtkImageData::SafeDownCast(input);
        identifiers.setDims(data->GetDimensions());
        identifiers.setSpacing(data->GetSpacing());
        status = identifiers.executeImageData();
        break;
      }
      default: {
        this->printErr("Unable to triangulate `"
                       + std::string(input->GetClassName()) + "`");
      }
    }
  } else {
    status = identifiers.executeSequential();
  }

  if(status < 1) {
    printErr("Global identifier generation failed");
    return -1;
  }

  // Add VTK objects to the data set

  input->GetPointData()->SetGlobalIds(vtkVertexIdentifiers);
  input->GetCellData()->SetGlobalIds(vtkCellIdentifiers);
  return 0;
}

void ttkAlgorithm::MPIGhostPipelinePreconditioning(vtkDataSet *input) {

  vtkNew<vtkGhostCellsGenerator> generator;
  if(ttk::isRunningWithMPI()
     && (!input->HasAnyGhostCells()
         && ((input->GetPointData()->GetArray("RankArray") == nullptr)
             || (input->GetCellData()->GetArray("RankArray") == nullptr)))) {
    generator->SetInputData(input);
    generator->BuildIfRequiredOff();
    generator->SetNumberOfGhostLayers(1);
    generator->Update();
    input->ShallowCopy(generator->GetOutputDataObject(0));
    input->GetPointData()->AddArray(
      generator->GetOutputDataObject(0)->GetGhostArray(0));
    input->GetCellData()->AddArray(
      generator->GetOutputDataObject(0)->GetGhostArray(1));
  }
}

void ttkAlgorithm::MPIPipelinePreconditioning(
  vtkDataSet *input,
  std::vector<int> &neighbors,
  ttk::Triangulation *triangulation) {

  ttk::SimplexId vertexNumber = input->GetNumberOfPoints();
  ttk::SimplexId cellNumber = input->GetNumberOfCells();
  if((input->GetDataObjectType() == VTK_POLY_DATA
      || input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)) {
    if((input->GetCellData()->GetGlobalIds() == nullptr)
       || (input->GetPointData()->GetGlobalIds() == nullptr)) {
      printWrn("Up to Paraview 5.10.1, bugs have been found in VTK for the "
               "distribution of Unstructured Grids and Poly Data");
      printWrn("As a consequence, the generation of Global ids is "
               "incorrect for those simplices");
      printWrn(
        "Beware when using TTK for such data set types, some results may "
        "be false");
    }
    if((input->GetPointData()->GetArray("RankArray") == nullptr)
       || (input->GetCellData()->GetArray("RankArray") == nullptr)) {
      printWrn("Up to Paraview 5.10.1, bugs have been found in VTK for the "
               "distribution of Unstructured Grids and Poly Data");
      printWrn("As a consequence, the generation of RankArray is "
               "incorrect for those simplices");
      printWrn(
        "Beware when using TTK for such data set types, some results may "
        "be false");
    }
  }

  // Get the neighbor ranks
  std::vector<int> &neighborRanks{
    triangulation != nullptr ? triangulation->getNeighborRanks() : neighbors};

  double *boundingBox = input->GetBounds();
  if(triangulation != nullptr) {
    triangulation->createMetaGrid(boundingBox);
  }

  if(neighborRanks.empty()) {
    ttk::preconditionNeighborsUsingBoundingBox(boundingBox, neighborRanks);
  }

  // Checks if global ids are valid
  ttk::LongSimplexId *globalPointIds = ttkUtils::GetPointer<ttk::LongSimplexId>(
    input->GetPointData()->GetGlobalIds());
  ttk::LongSimplexId *globalCellIds = ttkUtils::GetPointer<ttk::LongSimplexId>(
    input->GetCellData()->GetGlobalIds());

  bool pointValidity{false};
  bool cellValidity{false};
  if((triangulation != nullptr
      && (triangulation->getType() == ttk::Triangulation::Type::EXPLICIT
          || triangulation->getType() == ttk::Triangulation::Type::COMPACT))
     || triangulation == nullptr) {
    if(globalPointIds != nullptr) {
      unsigned char *ghostPoints = ttkUtils::GetPointer<unsigned char>(
        input->GetPointData()->GetArray("vtkGhostType"));
      int *vertexRankArray = ttkUtils::GetPointer<int>(
        input->GetPointData()->GetArray("RankArray"));
      pointValidity = checkGlobalIdValidity(
        globalPointIds, vertexNumber, ghostPoints, vertexRankArray);
    }
    if(pointValidity && globalCellIds != nullptr) {

      unsigned char *ghostCells = ttkUtils::GetPointer<unsigned char>(
        input->GetCellData()->GetArray("vtkGhostType"));
      int *cellRankArray = ttkUtils::GetPointer<int>(
        input->GetCellData()->GetArray("RankArray"));
      cellValidity = checkGlobalIdValidity(
        globalCellIds, cellNumber, ghostCells, cellRankArray);
    }
  } else {
    pointValidity = true;
    cellValidity = true;
  }

  // If the global ids are not valid, they are computed again
  if(!pointValidity || !cellValidity) {
    if(triangulation != nullptr) {
      if(triangulation->getType() == ttk::Triangulation::Type::EXPLICIT
         || triangulation->getType() == ttk::Triangulation::Type::COMPACT) {
        this->GenerateGlobalIds(
          input, triangulation->getVertexGlobalIdMap(), neighborRanks);
      }
    } else {
      std::unordered_map<ttk::SimplexId, ttk::SimplexId> vertGtoL{};
      this->GenerateGlobalIds(input, vertGtoL, neighborRanks);
    }
  }
}

void ttkAlgorithm::MPITriangulationPreconditioning(
  ttk::Triangulation *triangulation, vtkDataSet *input) {

  const auto pd{input->GetPointData()};
  if(pd == nullptr) {
    triangulation->printWrn("No point data on input object");
  } else {
    // provide "GlobalPointIds" & "vtkGhostType" point data arrays
    // to the triangulation
    triangulation->setVertsGlobalIds(
      ttkUtils::GetPointer<ttk::LongSimplexId>(pd->GetGlobalIds()));
    triangulation->setVertexGhostArray(
      ttkUtils::GetPointer<unsigned char>(pd->GetArray("vtkGhostType")));
    int *vertexRankArray = ttkUtils::GetPointer<int>(pd->GetArray("RankArray"));
    if(vertexRankArray != nullptr) {
      triangulation->setVertexRankArray(vertexRankArray);
    }
    triangulation->preconditionDistributedVertices();
  }

  const auto cd{input->GetCellData()};
  if(cd == nullptr) {
    triangulation->printWrn("No cell data on input object");
  } else {
    // provide "GlobalCellIds" & "vtkGhostType" cell data arrays to
    // the triangulation
    triangulation->setCellsGlobalIds(
      ttkUtils::GetPointer<ttk::LongSimplexId>(cd->GetGlobalIds()));
    triangulation->setCellGhostArray(
      ttkUtils::GetPointer<unsigned char>(cd->GetArray("vtkGhostType")));
    int *cellRankArray = ttkUtils::GetPointer<int>(cd->GetArray("RankArray"));
    if(cellRankArray != nullptr) {
      triangulation->setCellRankArray(cellRankArray);
    }
    triangulation->preconditionDistributedCells();
  }
}

#endif // TTK_ENABLE_MPI

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
#ifdef TTK_ENABLE_MPI
    if(ttk::hasInitializedMPI()) {
      if(this->updateMPICommunicator(vtkDataSet::GetData(inputVector[0], 0))) {
        return 1;
      };
    }
#endif // TTK_ENABLE_MPI
    return this->RequestData(request, inputVector, outputVector);
  }

  this->printErr("Unsupported pipeline pass:");
  request->Print(cout);

  return 0;
}
