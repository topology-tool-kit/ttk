#include <OrderDisambiguation.h>
#include <mpi.h>
#include <ttkArrayPreconditioning.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <regex>
#include <unordered_map>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkArrayPreconditioning);
using IT = long long int;

ttkArrayPreconditioning::ttkArrayPreconditioning() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->setDebugMsgPrefix("ArrayPreconditioning");

  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this, &ttkArrayPreconditioning::Modified);
}

int ttkArrayPreconditioning::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkArrayPreconditioning::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

void ttkArrayPreconditioning::ReceiveAndAddToVector(
  MPI_Datatype mpi_values,
  int rankFrom,
  int structTag,
  int intTag,
  std::vector<std::vector<ttk::value>> &unsortedReceivedValues) {
  std::vector<ttk::value> receivedValues;
  // be prepared to receive burstsize of elements, resize after receiving to the
  // correct size
  int amount;
  MPI_Recv(
    &amount, 1, MPI_INT, rankFrom, intTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  receivedValues.resize(amount, {0, 0, 0});
  // this->printMsg("Receiving " + std::to_string(amount) + " values from rank "
  // + std::to_string(rankFrom));
  MPI_Recv(receivedValues.data(), amount, mpi_values, rankFrom, structTag,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  unsortedReceivedValues[rankFrom] = receivedValues;
}

int ttkArrayPreconditioning::RequestData(vtkInformation *ttkNotUsed(request),
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
  int nVertices = input->GetNumberOfPoints();

  std::vector<vtkDataArray *> scalarArrays{};

  if(SelectFieldsWithRegexp) {
    // select all input point data arrays whose name is matching the regexp
    const auto n = pointData->GetNumberOfArrays();
    for(int i = 0; i < n; ++i) {
      auto array = pointData->GetArray(i);
      if(array != nullptr && array->GetName() != nullptr
         && std::regex_match(array->GetName(), std::regex(RegexpString))) {
        scalarArrays.emplace_back(array);
      }
    }
  } else {
    // get all selected input point data arrays
    for(int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
      auto array = pointData->GetArray(i);
      if(array != nullptr && array->GetName() != nullptr
         && ArraySelection->ArrayIsEnabled(array->GetName())) {
        scalarArrays.emplace_back(array);
      }
    }
  }

  auto vtkGlobalPointIds = pointData->GetGlobalIds();
  auto vtkGhostCells = pointData->GetArray("vtkGhostType");
  auto rankArray = pointData->GetArray("RankArray");
  if(vtkGlobalPointIds != nullptr && vtkGhostCells != nullptr
     && rankArray != nullptr) {
    MPI_Datatype mpi_values;
    const int nitems = 4;
    int blocklengths[4] = {1, 1, 1, 1};
    MPI_Datatype MIT = MPI_LONG_LONG_INT;
    MPI_Datatype types[4] = {MPI_FLOAT, MIT, MIT, MIT};
    MPI_Aint offsets[4];
    offsets[0] = offsetof(ttk::value, scalar);
    offsets[1] = offsetof(ttk::value, globalId);
    offsets[2] = offsetof(ttk::value, localId);
    offsets[3] = offsetof(ttk::value, ordering);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_values);
    MPI_Type_commit(&mpi_values);

    int numProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int intTag = 100;
    int structTag = 101;
    int boolTag = 102;
    if(rank == 0)
      this->printMsg("Global Point Ids, Ghost Cells and RankArray exist, "
                     "therefore we are in distributed mode!");
    this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank "
                   + std::to_string(rank));
    MPI_Barrier(MPI_COMM_WORLD);
    // add the order array for every scalar array, except the ghostcells and the
    // global ids
    for(auto scalarArray : scalarArrays) {
      std::string arrayName = std::string(scalarArray->GetName());
      if(arrayName != "GlobalPointIds" && arrayName != "vtkGhostType"
         && arrayName != "RankArray") {

        if(rank == 0)
          this->printMsg("Arrayname: " + arrayName);
        this->printMsg("#Points in Rank " + std::to_string(rank) + ": "
                       + std::to_string(nVertices));
        ttk::Timer fillAndSortTimer;
        std::vector<ttk::value> sortingValues;
        std::vector<IT> gidsToGetVector;
        std::unordered_map<IT, IT> gidToLidMap;
        tie(sortingValues, gidsToGetVector, gidToLidMap) = ttk::populateVector(
          nVertices, ttkUtils::GetPointer<float>(scalarArray),
          ttkUtils::GetPointer<IT>(vtkGlobalPointIds),
          ttkUtils::GetPointer<char>(vtkGhostCells));

        // sort the scalar array distributed first by the scalar value itself,
        // then by the global id
        ttk::sortVerticesDistributed(sortingValues);
        this->printMsg("#Unique Points in Rank " + std::to_string(rank) + ": "
                       + std::to_string(sortingValues.size()));
        this->printMsg("#Ghostpoints in Rank " + std::to_string(rank) + ": "
                       + std::to_string(gidsToGetVector.size()));
        // when all are done sorting, rank 0 requests the highest values and
        // merges them
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0) {
          this->printMsg(
            "Filling vector and sorting for each rank done, starting merge.", 1,
            fillAndSortTimer.getElapsedTime());
        }

        std::vector<IT> orderedValuesForRank;
        std::vector<ttk::value> finalValues;
        ttk::Timer mergeTimer;
        size_t totalSize = 0;
        if(rank == 0) {
          this->printMsg("Rank 0 starts merging");
          // get the nVertices from each rank, add them to get the complete size
          // of the dataset
          for(int i = 0; i < numProcs; i++) {
            if(i == 0) {
              totalSize += sortingValues.size();
            } else {
              IT receivedSize;
              MPI_Recv(&receivedSize, 1, MIT, i, intTag, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
              totalSize += receivedSize;
            }
          }
          this->printMsg("Total amount of distributed points: "
                         + std::to_string(totalSize));
          IT currentOrder = totalSize - 1;
          std::vector<std::vector<ttk::value>> unsortedReceivedValues;
          unsortedReceivedValues.resize(numProcs);
          std::vector<std::vector<IT>> orderResendValues;
          orderResendValues.resize(numProcs);

          // receive the first batch of values
          for(int i = 0; i < numProcs; i++) {
            if(i == 0) {
              std::vector<ttk::value> ownValues
                = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
              unsortedReceivedValues[i] = ownValues;
            } else {
              this->ReceiveAndAddToVector(
                mpi_values, i, structTag, intTag, unsortedReceivedValues);
            }
          }

          while(finalValues.size() < totalSize) {
            // take the current maximum scalar over all ranks
            int rankIdOfMaxScalar = -1;
            float maxScalar = -FLT_MAX;
            IT maxGId = -1;
            for(int i = 0; i < numProcs; i++) {
              if(unsortedReceivedValues[i].size() > 0) {
                int thisId = i;
                float thisScalar = unsortedReceivedValues[i].back().scalar;
                IT thisGId = unsortedReceivedValues[i].back().globalId;
                if(thisScalar > maxScalar
                   || (thisScalar == maxScalar && thisGId > maxGId)) {
                  maxScalar = thisScalar;
                  maxGId = thisGId;
                  rankIdOfMaxScalar = thisId;
                }
              }
            }
            if(rankIdOfMaxScalar == -1) {
              this->printMsg("FinalValues.size: "
                             + std::to_string(finalValues.size()));
              this->printMsg("All vectors are empty, but out final vector is "
                             "not complete yet. Either something went wrong or "
                             "some rank didn't send their values yet.");
              return 0;
            }

            // move the struct from the unsortedReceivedValues subvector to the
            // finalValues vector to get an ordering
            ttk::value currentValue
              = unsortedReceivedValues[rankIdOfMaxScalar].back();
            // we send the globalId and the the order via one send command,
            // therefore we need to check two concurrent values at the same time
            // later on
            // this->printMsg("Current GlobalId is " +
            // std::to_string(currentValue.globalId) + " from rank " +
            // std::to_string(rankIdOfMaxScalar));
            orderResendValues[rankIdOfMaxScalar].push_back(
              currentValue.globalId);
            orderResendValues[rankIdOfMaxScalar].push_back(currentOrder);
            currentOrder--;
            finalValues.push_back(currentValue);
            unsortedReceivedValues[rankIdOfMaxScalar].pop_back();
            if(unsortedReceivedValues[rankIdOfMaxScalar].size() == 0) {
              // this->printMsg("Vector for Rank " +
              // std::to_string(rankIdOfMaxScalar) + " is empty, we either need
              // to receive more or are done");

              if(rankIdOfMaxScalar == 0) {
                // append the ordered values to the correct vector
                orderedValuesForRank.insert(
                  orderedValuesForRank.end(),
                  orderResendValues[rankIdOfMaxScalar].begin(),
                  orderResendValues[rankIdOfMaxScalar].end());
                orderResendValues[rankIdOfMaxScalar].clear();
                if(sortingValues.size() > 0) {
                  std::vector<ttk::value> ownValues
                    = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
                  unsortedReceivedValues[rankIdOfMaxScalar] = ownValues;
                } else {
                  this->printMsg("We are done with rank 0!");
                }
              } else {
                // receive more values from rank, send ordering to the rank
                // send to the finished rank that we want more
                MPI_Send(orderResendValues[rankIdOfMaxScalar].data(),
                         orderResendValues[rankIdOfMaxScalar].size(), MIT,
                         rankIdOfMaxScalar, intTag, MPI_COMM_WORLD);
                orderResendValues[rankIdOfMaxScalar].clear();
                // check if there are more values to be received. If so, receive
                // them
                bool moreVals = false;
                MPI_Recv(&moreVals, 1, MPI_CXX_BOOL, rankIdOfMaxScalar, boolTag,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(moreVals) {
                  // this->printMsg("ANOTHER ONE " +
                  // std::to_string(rankIdOfMaxScalar));
                  this->ReceiveAndAddToVector(mpi_values, rankIdOfMaxScalar,
                                              structTag, intTag,
                                              unsortedReceivedValues);
                } else {
                  this->printMsg("We are done with rank "
                                 + std::to_string(rankIdOfMaxScalar));
                }
              }
            }
          }

          this->printMsg("Finished with sorting, max value is "
                         + std::to_string(finalValues[0].scalar)
                         + ", min value is "
                         + std::to_string(finalValues.back().scalar));
        } else {
          IT nValues = sortingValues.size();
          MPI_Send(&nValues, 1, MIT, 0, intTag, MPI_COMM_WORLD);

          // send the next burstsize values and then wait for an answer from the
          // root rank
          while(sortingValues.size() > 0) {
            std::vector<ttk::value> sendValues
              = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
            // this->printMsg("Rank " + std::to_string(rank) + " sending " +
            // std::to_string(sendValues.size()) + " values to rank 0, first
            // scalar is " + std::to_string(sendValues.data()[0].scalar));
            int size = sendValues.size();
            MPI_Send(&size, 1, MPI_INT, 0, intTag, MPI_COMM_WORLD);
            MPI_Send(sendValues.data(), size, mpi_values, 0, structTag,
                     MPI_COMM_WORLD);
            std::vector<IT> receivedValues;

            // be prepared to receive burstsize of elements, resize after
            // receiving to the correct size
            receivedValues.resize(BurstSize * 2);
            MPI_Status status;
            int amount;
            MPI_Recv(receivedValues.data(), BurstSize * 2, MIT, 0, intTag,
                     MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MIT, &amount);
            receivedValues.resize(amount);
            orderedValuesForRank.insert(orderedValuesForRank.end(),
                                        receivedValues.begin(),
                                        receivedValues.end());

            // afterwards send to root if there are still values to be sent
            bool moreVals = sortingValues.size() > 0;
            MPI_Send(&moreVals, 1, MPI_CXX_BOOL, 0, boolTag, MPI_COMM_WORLD);
          }
        }

        // all ranks do the following
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0) {
          this->printMsg(
            "Merging done and results sent to ranks, ranks are getting order "
            "for ghost cells and constructing order array.",
            1, mergeTimer.getElapsedTime());
        }
        this->printMsg("#Orders received for Rank " + std::to_string(rank)
                       + ": " + std::to_string(orderedValuesForRank.size()));

        /*ttk::Timer sendTimer;
        MPI_Bcast(finalValues.data(), totalSize, mpi_values, 0, MPI_COMM_WORLD);
        this->printMsg("Sent results, generating order array for rank " +
        std::to_string(rank), 1, sendTimer.getElapsedTime());
        */
        ttk::Timer orderTimer;
        vtkNew<ttkSimplexIdTypeArray> orderArray{};
        orderArray->SetName(this->GetOrderArrayName(scalarArray).data());
        orderArray->SetNumberOfComponents(1);
        orderArray->SetNumberOfTuples(nVertices);

        ttk::buildArrayForReceivedData(
          orderedValuesForRank.size(), orderedValuesForRank.data(), gidToLidMap,
          ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
          this->threadNumber_);

        // we still need to get the data for the ghostcells from their ranks
        auto neighbors = ttk::getNeighbors(
          nVertices, rank, ttkUtils::GetPointer<int>(rankArray));
        this->printMsg("Rank " + std::to_string(rank)
                       + " #Neighbors: " + std::to_string(neighbors.size()));
        auto gIdForNeighbors = ttk::getGIdForNeighbors(
          gidsToGetVector.size(), neighbors.size(), gidsToGetVector, neighbors,
          gidToLidMap, ttkUtils::GetPointer<int>(rankArray));
        MPI_Request req;
        std::vector<size_t> sizesToSend;
        sizesToSend.resize(gIdForNeighbors.size());
        for(size_t i = 0; i < gIdForNeighbors.size(); i++) {
          auto iter = neighbors.begin();
          std::advance(iter, i);
          int rankToSend = *iter;
          std::vector<IT> gIdToSend = gIdForNeighbors[i];
          sizesToSend[i] = gIdToSend.size();
          // we need to send to all neighbors how many we will send them,
          // and then send the values to neighbors where the amount > 0
          // first we Isend everything we need to the neighbors
          MPI_Isend(
            &sizesToSend[i], 1, MIT, rankToSend, intTag, MPI_COMM_WORLD, &req);
          MPI_Request_free(&req);
          if(sizesToSend[i] > 0) {
            MPI_Isend(gIdForNeighbors[i].data(), sizesToSend[i], MIT,
                      rankToSend, intTag, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
          }
        }

        // then we blockingly receive from our neighbors
        std::vector<std::vector<IT>> ordersToSend;
        ordersToSend.resize(neighbors.size());
        for(size_t i = 0; i < neighbors.size(); i++) {
          auto iter = neighbors.begin();
          std::advance(iter, i);
          int rankToRecv = *iter;
          size_t sizeToRecv;
          MPI_Recv(&sizeToRecv, 1, MIT, rankToRecv, intTag, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);

          if(sizeToRecv > 0) {
            std::vector<IT> gIdToRecv;
            gIdToRecv.resize(sizeToRecv);
            MPI_Recv(gIdToRecv.data(), sizeToRecv, MIT, rankToRecv, intTag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // prepare what they need
            ordersToSend[i] = ttk::getOrderForGIds(
              gIdToRecv.size(), gIdToRecv.data(), gidToLidMap,
              ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
              this->threadNumber_);

            // Isend it to them
            MPI_Isend(ordersToSend[i].data(), ordersToSend[i].size(), MIT,
                      rankToRecv, intTag, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
          }
        }

        // and blockingly receive the order from our neighbors
        for(size_t i = 0; i < gIdForNeighbors.size(); i++) {
          auto iter = neighbors.begin();
          std::advance(iter, i);
          int rankToRecvOrder = *iter;
          auto sizeToRecvOrder = gIdForNeighbors[i].size() * 2;
          std::vector<IT> ordersToRecv;
          ordersToRecv.resize(sizeToRecvOrder);
          if(sizeToRecvOrder > 0) {
            MPI_Recv(ordersToRecv.data(), sizeToRecvOrder, MIT, rankToRecvOrder,
                     intTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // and add it to our orderarray
            ttk::buildArrayForReceivedData(
              ordersToRecv.size(), ordersToRecv.data(), gidToLidMap,
              ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
              this->threadNumber_);
          }
        }

        // this->printMsg("Neighbors for rank " + std::to_string(rank) + ": " +
        // std::to_string(neighbors));
        /*for (int i = 0; i < numProcs; i++){
          MPI_Send()
        }*/

        output->GetPointData()->AddArray(orderArray);
        this->printMsg("Generated order array for scalar array `"
                         + std::string{scalarArray->GetName()} + "', rank "
                         + std::to_string(rank),
                       1, orderTimer.getElapsedTime());

        /*
        std::unordered_map<float, int> orderMap =
        ttk::buildOrderMap(finalValues, totalSize);
        // every rank now has an orderedValuesForRank array with the points
        sorted in descending order and their correct order
        // now we need to transform this to a correct vtk orderarray and append
        it



        ttk::buildOrderArray(nVertices,
                            ttkUtils::GetPointer<float>(scalarArray),
                            orderMap,
                            ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
                            this->threadNumber_);


        */
      }
    }
    this->printMsg("Preconditioned selected scalar arrays", 1.0,
                   tm.getElapsedTime(), this->threadNumber_);
    return 1;
  }

  for(auto scalarArray : scalarArrays) {
    vtkNew<ttkSimplexIdTypeArray> orderArray{};
    orderArray->SetName(this->GetOrderArrayName(scalarArray).data());
    orderArray->SetNumberOfComponents(1);
    orderArray->SetNumberOfTuples(nVertices);

    switch(scalarArray->GetDataType()) {
      vtkTemplateMacro(ttk::sortVertices(
        nVertices, static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)),
        static_cast<int *>(nullptr),
        static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(orderArray)),
        this->threadNumber_));
    }

    output->GetPointData()->AddArray(orderArray);
    this->printMsg("Generated order array for scalar array `"
                   + std::string{scalarArray->GetName()} + "'");
  }

  this->printMsg("Preconditioned selected scalar arrays", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
