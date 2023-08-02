#include <string>
#include <ttkIntegralLines.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <ArrayLinkedList.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

#include <array>

vtkStandardNewMacro(ttkIntegralLines);

ttkIntegralLines::ttkIntegralLines() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkIntegralLines::~ttkIntegralLines() = default;

int ttkIntegralLines::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int ttkIntegralLines::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");

  return 1;
}

template <typename triangulationType>
int ttkIntegralLines::getTrajectories(
  vtkDataSet *input,
  const triangulationType *triangulation,
  const std::vector<
    ttk::ArrayLinkedList<ttk::intgl::IntegralLine, INTEGRAL_LINE_TABULAR_SIZE>>
    &integralLines,
#ifdef TTK_ENABLE_MPI
  const std::vector<ttk::SimplexId> &globalVertexId,
  const std::vector<ttk::SimplexId> &globalCellId,
#endif
  vtkUnstructuredGrid *output) {
  if(input == nullptr || output == nullptr
     || input->GetPointData() == nullptr) {
    this->printErr("Null pointers in getTrajectories parameters");
    return 0;
  }

  vtkNew<vtkUnstructuredGrid> ug{};
  vtkNew<vtkPoints> pts{};
  vtkNew<vtkDoubleArray> dist{};
  vtkNew<vtkIdTypeArray> identifier{};
#ifdef TTK_ENABLE_MPI
  vtkNew<vtkIdTypeArray> vtkEdgeIdentifiers{};
  vtkNew<vtkIdTypeArray> vtkVertexGlobalIdArray{};
  vtkNew<vtkIntArray> vtkVertexRankArray{};
  vtkNew<vtkIntArray> vtkEdgeRankArray{};
#endif
  vtkNew<vtkIdTypeArray> vtkForkIdentifiers{};
  vtkNew<vtkUnsignedCharArray> outputMaskField{};

  outputMaskField->SetNumberOfComponents(1);
  outputMaskField->SetName(ttk::MaskScalarFieldName);

  dist->SetNumberOfComponents(1);
  dist->SetName("DistanceFromSeed");
  identifier->SetNumberOfComponents(1);
  identifier->SetName("SeedIdentifier");
  vtkForkIdentifiers->SetNumberOfComponents(1);
  vtkForkIdentifiers->SetName("ForkIdentifiers");

#ifdef TTK_ENABLE_MPI
  vtkVertexGlobalIdArray->SetNumberOfComponents(1);
  vtkVertexGlobalIdArray->SetName("GlobalPointIds");
  vtkEdgeIdentifiers->SetNumberOfComponents(1);
  vtkEdgeIdentifiers->SetName("GlobalCellIds");
  vtkEdgeRankArray->SetNumberOfComponents(1);
  vtkEdgeRankArray->SetName("RankArray");
  vtkVertexRankArray->SetNumberOfComponents(1);
  vtkVertexRankArray->SetName("RankArray");
#endif
  const auto numberOfArrays = input->GetPointData()->GetNumberOfArrays();

  std::vector<vtkDataArray *> scalarArrays{};
  scalarArrays.reserve(numberOfArrays);
  for(int k = 0; k < numberOfArrays; ++k) {
    const auto a = input->GetPointData()->GetArray(k);
    if(a->GetNumberOfComponents() == 1) {
      // only keep scalar arrays
      scalarArrays.push_back(a);
    }
  }

  std::vector<vtkSmartPointer<vtkDataArray>> inputScalars(scalarArrays.size());
  for(size_t k = 0; k < scalarArrays.size(); ++k) {
    inputScalars[k]
      = vtkSmartPointer<vtkDataArray>::Take(scalarArrays[k]->NewInstance());
    inputScalars[k]->SetNumberOfComponents(1);
    inputScalars[k]->SetName(scalarArrays[k]->GetName());
  }
  std::array<float, 3> p;
  std::array<vtkIdType, 2> ids;
#ifdef TTK_ENABLE_MPI
  ttk::SimplexId vertexCounter = 0;
  ttk::SimplexId edgeCounter = 0;
#endif
  for(int thread = 0; thread < threadNumber_; thread++) {
    auto integralLine = integralLines[thread].list_.begin();
    while(integralLine != integralLines[thread].list_.end()) {
      for(int i = 0; i < INTEGRAL_LINE_TABULAR_SIZE; i++) {
        if(integralLine->at(i).trajectory.size() > 0) {
          ttk::SimplexId vertex = integralLine->at(i).trajectory.at(0);
          triangulation->getVertexPoint(vertex, p[0], p[1], p[2]);
          ids[0] = pts->InsertNextPoint(p.data());
          // distanceScalars
          dist->InsertNextTuple1(integralLine->at(i).distanceFromSeed.at(0));
#ifdef TTK_ENABLE_MPI
          outputMaskField->InsertNextTuple1(0);
          vtkVertexGlobalIdArray->InsertNextTuple1(
            globalVertexId.at(vertexCounter));
          vertexCounter++;
          vtkVertexRankArray->InsertNextTuple1(
            triangulation->getVertexRank(vertex));
#else
          outputMaskField->InsertNextTuple1(0);
#endif
          identifier->InsertNextTuple1(integralLine->at(i).seedIdentifier);
          vtkForkIdentifiers->InsertNextTuple1(
            integralLine->at(i).forkIdentifier);
          // inputScalars
          for(size_t k = 0; k < scalarArrays.size(); ++k) {
            inputScalars[k]->InsertNextTuple1(
              scalarArrays[k]->GetTuple1(vertex));
          }
          for(size_t j = 1; j < integralLine->at(i).trajectory.size(); ++j) {
            vertex = integralLine->at(i).trajectory.at(j);
#ifdef TTK_ENABLE_MPI
            vtkVertexGlobalIdArray->InsertNextTuple1(
              globalVertexId.at(vertexCounter));
            vertexCounter++;
            vtkEdgeIdentifiers->InsertNextTuple1(globalCellId.at(edgeCounter));
            edgeCounter++;
            vtkEdgeRankArray->InsertNextTuple1(triangulation->getVertexRank(
              integralLine->at(i).trajectory.at(j - 1)));
            vtkVertexRankArray->InsertNextTuple1(
              triangulation->getVertexRank(vertex));
#endif
            outputMaskField->InsertNextTuple1(1);
            vtkForkIdentifiers->InsertNextTuple1(
              integralLine->at(i).forkIdentifier);
            triangulation->getVertexPoint(vertex, p[0], p[1], p[2]);
            ids[1] = pts->InsertNextPoint(p.data());
            // distanceScalars
            dist->InsertNextTuple1(integralLine->at(i).distanceFromSeed.at(j));
            identifier->InsertNextTuple1(integralLine->at(i).seedIdentifier);
            // inputScalars
            for(unsigned int k = 0; k < scalarArrays.size(); ++k)
              inputScalars[k]->InsertNextTuple1(
                scalarArrays[k]->GetTuple1(vertex));
            ug->InsertNextCell(VTK_LINE, 2, ids.data());
            // iteration
            ids[0] = ids[1];
          }
          outputMaskField->SetTuple1(
            outputMaskField->GetNumberOfTuples() - 1, 0);
        } else {
          break;
        }
      }
      integralLine++;
    }
  }

  ug->SetPoints(pts);
  ug->GetPointData()->AddArray(dist);
  ug->GetPointData()->AddArray(identifier);
  ug->GetPointData()->AddArray(vtkForkIdentifiers);
  ug->GetPointData()->AddArray(outputMaskField);

  for(unsigned int k = 0; k < scalarArrays.size(); ++k) {
    ug->GetPointData()->AddArray(inputScalars[k]);
  }
#ifdef TTK_ENABLE_MPI
  ug->GetPointData()->AddArray(vtkVertexRankArray);
  ug->GetCellData()->AddArray(vtkEdgeRankArray);
  ug->GetCellData()->SetGlobalIds(vtkEdgeIdentifiers);
  ug->GetPointData()->SetGlobalIds(vtkVertexGlobalIdArray);
#endif
  output->ShallowCopy(ug);

  return 1;
}

int ttkIntegralLines::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  vtkDataSet *domain = vtkDataSet::GetData(inputVector[0], 0);
  vtkPointSet *seeds = vtkPointSet::GetData(inputVector[1], 0);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(domain);
  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, domain);

  int const keepGoing = checkEmptyMPIInput<ttk::Triangulation>(triangulation);
  if(keepGoing < 2) {
    return keepGoing;
  }
  vtkDataArray *inputOffsets
    = this->GetOrderArray(domain, 0, 1, ForceInputOffsetScalarField);

  const ttk::SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();
  this->setVertexNumber(numberOfPointsInDomain);
  int numberOfPointsInSeeds = seeds->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  int totalSeeds;
#else
#ifdef TTK_ENABLE_MPI
  int totalSeeds;
#endif
#endif

#ifdef TTK_ENABLE_MPI_TIME
  ttk::Timer t_mpi;
  ttk::startMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
#endif
#ifdef TTK_ENABLE_MPI
  // Necessary when using MPI
  std::vector<ttk::SimplexId> inputIdentifiers{};
  if(ttk::MPIsize_ > 1) {
    MPI_Reduce(&numberOfPointsInSeeds, &totalSeeds, 1, MPI_INTEGER, MPI_SUM, 0,
               ttk::MPIcomm_);
    int isDistributed;

    if(ttk::MPIrank_ == 0) {
      isDistributed = numberOfPointsInSeeds != totalSeeds;
    }
    MPI_Bcast(&isDistributed, 1, MPI_INTEGER, 0, ttk::MPIcomm_);
    MPI_Bcast(&totalSeeds, 1, MPI_INTEGER, 0, ttk::MPIcomm_);

    if(!isDistributed) {
      this->setGlobalElementCounter(totalSeeds);
      vtkDataArray *globalSeedsId;
      if(ttk::MPIrank_ == 0) {
        globalSeedsId = seeds->GetPointData()->GetArray("GlobalPointIds");
      } else {
        globalSeedsId = vtkDataArray::CreateDataArray(VTK_ID_TYPE);
      }

      if(ttk::MPIrank_ != 0) {
        globalSeedsId->SetNumberOfComponents(1);
        globalSeedsId->SetNumberOfTuples(totalSeeds);
      }
      ttk::LongSimplexId id = 0;
      MPI_Bcast(ttkUtils::GetPointer<ttk::LongSimplexId>(globalSeedsId),
                totalSeeds, ttk::getMPIType(id), 0, ttk::MPIcomm_);
      ttk::SimplexId localId = -1;
      for(int i = 0; i < totalSeeds; i++) {
        localId = triangulation->getVertexLocalId(globalSeedsId->GetTuple1(i));

        if(localId != -1
           && triangulation->getVertexRank(localId) == ttk::MPIrank_) {
          inputIdentifiers.push_back(localId);
        }
      }
      numberOfPointsInSeeds = inputIdentifiers.size();
    } else {
      std::vector<ttk::SimplexId> idSpareStorage{};
      ttk::SimplexId *inputIdentifierGlobalId;
      inputIdentifierGlobalId = this->GetIdentifierArrayPtr(
        ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
        idSpareStorage);
      ttk::SimplexId localId = 0;
      for(int i = 0; i < numberOfPointsInSeeds; i++) {
        localId = triangulation->getVertexLocalId(inputIdentifierGlobalId[i]);
        if(localId != -1) {
          inputIdentifiers.push_back(localId);
        }
      }
      numberOfPointsInSeeds = inputIdentifiers.size();
      MPI_Allreduce(&numberOfPointsInSeeds, &totalSeeds, 1, MPI_INTEGER,
                    MPI_SUM, ttk::MPIcomm_);
      this->setGlobalElementCounter(totalSeeds);
    }
  } else {
    this->setGlobalElementCounter(numberOfPointsInSeeds);
    inputIdentifiers.resize(numberOfPointsInSeeds);
    totalSeeds = numberOfPointsInSeeds;
    std::vector<ttk::SimplexId> idSpareStorage{};
    ttk::SimplexId *inputIdentifierGlobalId;
    inputIdentifierGlobalId = this->GetIdentifierArrayPtr(
      ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
      idSpareStorage);
    for(int i = 0; i < numberOfPointsInSeeds; i++) {
      inputIdentifiers.at(i)
        = triangulation->getVertexLocalId(inputIdentifierGlobalId[i]);
    }
  }
#else
  std::vector<ttk::SimplexId> idSpareStorage{};
  ttk::SimplexId *identifiers = this->GetIdentifierArrayPtr(
    ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
    idSpareStorage);
  std::unordered_set<ttk::SimplexId> isSeed;
  for(ttk::SimplexId k = 0; k < numberOfPointsInSeeds; ++k) {
    isSeed.insert(identifiers[k]);
  }
  std::vector<ttk::SimplexId> inputIdentifiers(isSeed.begin(), isSeed.end());
#ifndef TTK_ENABLE_KAMIKAZE
  totalSeeds = inputIdentifiers.size();
#endif
  isSeed.clear();
#endif

  std::vector<
    ttk::ArrayLinkedList<ttk::intgl::IntegralLine, INTEGRAL_LINE_TABULAR_SIZE>>
    integralLines(
      threadNumber_, ttk::ArrayLinkedList<ttk::intgl::IntegralLine,
                                          INTEGRAL_LINE_TABULAR_SIZE>());

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSeedNumber(numberOfPointsInSeeds);
  this->setDirection(Direction);
  this->setInputScalarField(inputScalars->GetVoidPointer(0));
  this->setInputOffsets(ttkUtils::GetPointer<ttk::SimplexId>(inputOffsets));
  this->setVertexIdentifierScalarField(&inputIdentifiers);
  this->setOutputIntegralLines(&integralLines);
  this->preconditionTriangulation(triangulation);
  this->setChunkSize(
    std::max(std::max(std::min(1000, (int)numberOfPointsInSeeds),
                      (int)numberOfPointsInSeeds / (threadNumber_ * 100)),
             1));
#ifdef TTK_ENABLE_MPI
  std::vector<std::vector<std::vector<ttk::intgl::ElementToBeSent>>> toSend(
    ttk::MPIsize_);
  this->setNeighbors(triangulation->getNeighborRanks());
  if(ttk::MPIsize_ > 1) {
    toSend.resize(this->neighborNumber_);
    for(int i = 0; i < this->neighborNumber_; i++) {
      toSend[i].resize(this->threadNumber_);
      for(int j = 0; j < this->threadNumber_; j++) {
        toSend[i][j].reserve((int)numberOfPointsInSeeds * 0.005
                             / this->threadNumber_);
      }
    }
  }
  this->setToSend(&toSend);
  this->createMessageType();
#endif
#ifdef TTK_ENABLE_MPI_TIME
  double elapsedTime = ttk::endMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
  if(ttk::MPIrank_ == 0) {
    printMsg("Preparation performed using " + std::to_string(ttk::MPIsize_)
             + " MPI processes lasted :" + std::to_string(elapsedTime));
  }
#endif
#ifndef TTK_ENABLE_KAMIKAZE
  // field problem
  if(!inputScalars) {
    this->printErr("wrong scalar field.");
    return -1;
  }
  // field problem
  if(inputOffsets->GetDataType() != VTK_INT
     and inputOffsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("input offset field type not supported.");
    return -1;
  }

  // no points.
  if(numberOfPointsInDomain <= 0) {
    this->printErr("domain has no points.");
    return -1;
  }
  // no points.
  if(totalSeeds <= 0) {
    this->printErr("seeds have no points.");
    return -1;
  }
#endif
#ifdef TTK_ENABLE_MPI_TIME
  ttk::startMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
#endif
  int status = 0;
  ttkTemplateMacro(triangulation->getType(),
                   (status = this->execute<TTK_TT>(
                      static_cast<TTK_TT *>(triangulation->getData()))));

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(status != 0) {
    this->printErr("IntegralLines.execute() error code : "
                   + std::to_string(status));
    return 0;
  }
#endif
#ifdef TTK_ENABLE_MPI
  std::vector<ttk::SimplexId> globalVertexId;
  std::vector<ttk::SimplexId> globalCellId;
  ttkTemplateMacro(triangulation->getType(),
                   (getGlobalIdentifiers<TTK_TT>(
                     globalVertexId, globalCellId, integralLines,
                     static_cast<TTK_TT *>(triangulation->getData()))));
#ifdef TTK_ENABLE_MPI_TIME
  elapsedTime = ttk::endMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
  if(ttk::MPIrank_ == 0) {
    printMsg("Computation performed using " + std::to_string(ttk::MPIsize_)
             + " MPI processes lasted :" + std::to_string(elapsedTime));
  }
#endif
#endif
  // make the vtk trajectories
#ifdef TTK_ENABLE_MPI
  ttkTemplateMacro(triangulation->getType(),
                   (getTrajectories<TTK_TT>(
                     domain, static_cast<TTK_TT *>(triangulation->getData()),
                     integralLines, globalVertexId, globalCellId, output)));
#else
  ttkTemplateMacro(triangulation->getType(),
                   (getTrajectories<TTK_TT>(
                     domain, static_cast<TTK_TT *>(triangulation->getData()),
                     integralLines, output)));
#endif

  return (int)(status == 0);
}
