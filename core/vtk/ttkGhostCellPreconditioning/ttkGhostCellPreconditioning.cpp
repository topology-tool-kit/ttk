#include <ttkGhostCellPreconditioning.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <unordered_map>
#include <unordered_set>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkGhostCellPreconditioning);

ttkGhostCellPreconditioning::ttkGhostCellPreconditioning() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->setDebugMsgPrefix("GhostCellPreconditioning");

  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this, &ttkGhostCellPreconditioning::Modified);
}

int ttkGhostCellPreconditioning::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkGhostCellPreconditioning::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

// returns true if bounding boxes intersect, false if not
bool checkForIntersection(double *myBB, double *theirBB) {
  return !(
    myBB[0] > theirBB[1] // my left side is right of their right side
    || myBB[1] < theirBB[0] // my right side is left of their left side
    || myBB[2] > theirBB[3] // my bottom side is above their top side
    || myBB[3] < theirBB[2] // my top side is under their bottom side
    || myBB[4] > theirBB[5] // my front side is behind their back side
    || myBB[5] < theirBB[4] // my back side is in front of their front side
  );
}

int ttkGhostCellPreconditioning::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);
  ttk::Timer tm{};

  if(input == nullptr || output == nullptr) {
    return 0;
  }

  output->ShallowCopy(input);

  auto pointData = input->GetPointData();
  ttk::SimplexId nVertices = input->GetNumberOfPoints();
  this->printMsg("#Points: " + std::to_string(nVertices));

  auto vtkGlobalPointIds = pointData->GetGlobalIds();
  auto vtkGhostCells = pointData->GetArray("vtkGhostType");
  if(vtkGlobalPointIds != nullptr && vtkGhostCells != nullptr) {
#ifdef TTK_ENABLE_MPI
    if(ttk::isRunningWithMPI()) {
      if(ttk::MPIrank_ == 0)
        this->printMsg(
          "Global Point Ids and Ghost Cells exist, therefore we can continue!");
      this->printMsg("#Ranks " + std::to_string(ttk::MPIsize_)
                     + ", this is rank " + std::to_string(ttk::MPIrank_));
      double *boundingBox = input->GetBounds();
      std::vector<double *> rankBoundingBoxes(ttk::MPIsize_);
      rankBoundingBoxes[ttk::MPIrank_] = boundingBox;
      for(int r = 0; r < ttk::MPIsize_; r++) {
        if(r != ttk::MPIrank_)
          rankBoundingBoxes[r] = (double *)malloc(6 * sizeof(double));
        MPI_Bcast(rankBoundingBoxes[r], 6, MPI_DOUBLE, r, ttk::MPIcomm_);
      }

      double epsilon = 0.00001;
      // inflate our own bounding box by epsilon
      for(int i = 0; i < 6; i++) {
        if(i % 2 == 0)
          boundingBox[i] -= epsilon;
        if(i % 2 == 1)
          boundingBox[i] += epsilon;
      }
      std::vector<int> neighbors;
      for(int i = 0; i < ttk::MPIsize_; i++) {
        if(i != ttk::MPIrank_) {
          double *theirBoundingBox = rankBoundingBoxes[i];
          if(checkForIntersection(boundingBox, theirBoundingBox)) {
            neighbors.push_back(i);
          }
        }
      }

      MPI_Datatype MIT = ttk::getMPIType(static_cast<ttk::SimplexId>(0));
      vtkNew<vtkIntArray> rankArray{};
      rankArray->SetName("RankArray");
      rankArray->SetNumberOfComponents(1);
      rankArray->SetNumberOfTuples(nVertices);
      std::vector<ttk::SimplexId> currentRankUnknownIds;
      std::unordered_set<ttk::SimplexId> gIdSet;
      std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToLocalMap;

      for(int i = 0; i < nVertices; i++) {
        int ghostCellVal = vtkGhostCells->GetComponent(i, 0);
        ttk::SimplexId globalId = vtkGlobalPointIds->GetComponent(i, 0);
        if(ghostCellVal == 0) {
          // if the ghost cell value is 0, then this vertex mainly belongs to
          // this rank
          rankArray->SetComponent(i, 0, ttk::MPIrank_);
          gIdSet.insert(globalId);
        } else {
          // otherwise the vertex belongs to another rank and we need to find
          // out to which one this needs to be done by broadcasting the global
          // id and hoping for some other rank to answer
          currentRankUnknownIds.push_back(globalId);
          gIdToLocalMap[globalId] = i;
        }
      }

      ttk::SimplexId sizeOfCurrentRank = currentRankUnknownIds.size();
      std::vector<ttk::SimplexId> gIdsToSend;
      std::vector<ttk::SimplexId> receivedGlobals;
      receivedGlobals.resize(sizeOfCurrentRank);
      ttk::SimplexId sizeOfNeighbor;
      std::vector<ttk::SimplexId> neighborUnknownIds;
      for(int neighbor : neighbors) {
        // we first send the size and then all needed ids to the neighbor
        MPI_Sendrecv(&sizeOfCurrentRank, 1, MIT, neighbor, ttk::MPIrank_,
                     &sizeOfNeighbor, 1, MIT, neighbor, neighbor, ttk::MPIcomm_,
                     MPI_STATUS_IGNORE);
        neighborUnknownIds.resize(sizeOfNeighbor);
        gIdsToSend.reserve(sizeOfNeighbor);

        MPI_Sendrecv(currentRankUnknownIds.data(), sizeOfCurrentRank, MIT,
                     neighbor, ttk::MPIrank_, neighborUnknownIds.data(),
                     sizeOfNeighbor, MIT, neighbor, neighbor, ttk::MPIcomm_,
                     MPI_STATUS_IGNORE);

        // then we check if the needed globalid values are present in the local
        // globalid set if so, we send the rank value to the requesting rank
        for(ttk::SimplexId gId : neighborUnknownIds) {
          if(gIdSet.count(gId)) {
            // add the value to the vector which will be sent
            gIdsToSend.push_back(gId);
          }
        }
        MPI_Status status;
        int amount;

        MPI_Sendrecv(gIdsToSend.data(), gIdsToSend.size(), MIT, neighbor,
                     ttk::MPIrank_, receivedGlobals.data(),
                     currentRankUnknownIds.size(), MIT, neighbor, neighbor,
                     ttk::MPIcomm_, &status);

        MPI_Get_count(&status, MIT, &amount);
        receivedGlobals.resize(amount);

        for(ttk::SimplexId receivedGlobal : receivedGlobals) {
          ttk::SimplexId localVal = gIdToLocalMap[receivedGlobal];
          rankArray->SetComponent(localVal, 0, neighbor);
        }
        // cleanup
        gIdsToSend.clear();
        receivedGlobals.resize(sizeOfCurrentRank);
      }

      // free the communicator once we are done with everything MPI
      output->GetPointData()->AddArray(rankArray);

      this->printMsg("Preprocessed RankArray", 1.0, tm.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    } else {
      this->printMsg("Necessary arrays are present,  TTK is built with MPI "
                     "support, but not run with mpirun. Running sequentially.");
      return 0;
    }
#else
    this->printMsg(
      "Necessary arrays are present, but TTK is not built with MPI support");
    return 0;

#endif
  } else {
    this->printMsg("Necessary arrays are not present.");
    return 0;
  }
}
