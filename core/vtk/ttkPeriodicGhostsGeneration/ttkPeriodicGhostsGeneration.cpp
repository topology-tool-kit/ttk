#include <ttkPeriodicGhostsGeneration.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPeriodicGhostsGeneration);

ttkPeriodicGhostsGeneration::ttkPeriodicGhostsGeneration() {
  this->setDebugMsgPrefix("PeriodicGhostsGeneration");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}

int ttkPeriodicGhostsGeneration::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkPeriodicGhostsGeneration::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }
  return 0;
}

#ifdef TTK_ENABLE_MPI

int ttkPeriodicGhostsGeneration::RequestUpdateExtent(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *ttkNotUsed(outputVector)) {
  this->ComputeOutputExtent();
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  std::array<int, 6> outExtentWithoutGhost{outExtent_};
  for(int i = 0; i < 3; i++) {
    if(std::abs(outExtent_[2 * i] - outExtent_[2 * i + 1]) > spacing_[i] / 10) {
      if(!(localGlobalBounds_[2 * i].isBound == 1
           && localGlobalBounds_[2 * i].isBound
                == localGlobalBounds_[2 * i + 1].isBound)) {
        outExtentWithoutGhost[2 * i] += 1;
        outExtentWithoutGhost[2 * i + 1] -= 1;
      }
    }
  }
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
              outExtentWithoutGhost[0], outExtentWithoutGhost[1],
              outExtentWithoutGhost[2], outExtentWithoutGhost[3],
              outExtentWithoutGhost[4], outExtentWithoutGhost[5]);
  return 1;
}

int ttkPeriodicGhostsGeneration::RequestInformation(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVectors,
  vtkInformationVector *outputVector) {
  vtkInformation *inInfo = inputVectors[0]->GetInformationObject(0);
  inInfo->Get(
    vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExtent_.data());
  inInfo->Get(vtkDataObject::SPACING(), this->spacing_.data());
  inInfo->Get(vtkDataObject::ORIGIN(), this->origin_.data());
  this->ComputeOutputExtent();
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outExtent_.data(), 6);
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
  return 1;
}

int ttkPeriodicGhostsGeneration::ComputeOutputExtent() {

  // Retrieve the local extent for the current process
  vtkNew<vtkExtentTranslator> translator;
  translator->SetWholeExtent(inExtent_.data());
  translator->SetNumberOfPieces(ttk::MPIsize_);
  translator->SetPiece(ttk::MPIrank_);
  translator->PieceToExtent();
  int localExtent[6];
  translator->GetExtent(localExtent);

  // Compute the bounds of the local process
  double bounds[6] = {
    origin_[0] + localExtent[0] * spacing_[0],
    origin_[0] + localExtent[1] * spacing_[0],
    origin_[1] + localExtent[2] * spacing_[1],
    origin_[1] + localExtent[3] * spacing_[1],
    origin_[2] + localExtent[4] * spacing_[2],
    origin_[2] + localExtent[5] * spacing_[2],
  };
  // If the extent has never been computed or has changed
  if(!isOutputExtentComputed_
     || std::abs(boundsWithoutGhosts_[0] - bounds[0]) > spacing_[0] / 10
     || std::abs(boundsWithoutGhosts_[1] - bounds[1]) > spacing_[0] / 10
     || std::abs(boundsWithoutGhosts_[2] - bounds[2]) > spacing_[1] / 10
     || std::abs(boundsWithoutGhosts_[3] - bounds[3]) > spacing_[1] / 10
     || std::abs(boundsWithoutGhosts_[4] - bounds[4]) > spacing_[2] / 10
     || std::abs(boundsWithoutGhosts_[5] - bounds[5]) > spacing_[2] / 10) {
    boundsWithoutGhosts_
      = {bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]};
    // re-order tempGlobalBounds
    globalBounds_ = {origin_[0] + inExtent_[0] * spacing_[0],
                     origin_[0] + inExtent_[1] * spacing_[0],
                     origin_[1] + inExtent_[2] * spacing_[1],
                     origin_[1] + inExtent_[3] * spacing_[1],
                     origin_[2] + inExtent_[4] * spacing_[2],
                     origin_[2] + inExtent_[5] * spacing_[2]};

    ttk::preconditionNeighborsUsingBoundingBox(bounds, neighbors_);

    // Compute the 1D boundary matches of local process: if the process is on
    // the global boundary, then the center of its other dimension is computed

    for(int i = 0; i < 2; i++) {
      if(std::abs(globalBounds_[i] - boundsWithoutGhosts_[i])
         < spacing_[0] / 10) {
        localGlobalBounds_[i].isBound = 1;
        localGlobalBounds_[i].x = boundsWithoutGhosts_[i];
        localGlobalBounds_[i].y
          = (boundsWithoutGhosts_[2] + boundsWithoutGhosts_[3]) / 2;
        localGlobalBounds_[i].z
          = (boundsWithoutGhosts_[4] + boundsWithoutGhosts_[5]) / 2;
      }
    }

    for(int i = 0; i < 2; i++) {
      if(std::abs(globalBounds_[2 + i] - boundsWithoutGhosts_[2 + i])
         < spacing_[1] / 10) {
        localGlobalBounds_[2 + i].isBound = 1;
        localGlobalBounds_[2 + i].x
          = (boundsWithoutGhosts_[0] + boundsWithoutGhosts_[1]) / 2;
        localGlobalBounds_[2 + i].y = boundsWithoutGhosts_[2 + i];
        localGlobalBounds_[2 + i].z
          = (boundsWithoutGhosts_[4] + boundsWithoutGhosts_[5]) / 2;
      }
    }

    for(int i = 0; i < 2; i++) {
      if(std::abs(globalBounds_[4 + i] - boundsWithoutGhosts_[4 + i])
         < spacing_[2] / 10) {
        localGlobalBounds_[4 + i].isBound = 1;
        localGlobalBounds_[4 + i].x
          = (boundsWithoutGhosts_[0] + boundsWithoutGhosts_[1]) / 2;
        localGlobalBounds_[4 + i].y
          = (boundsWithoutGhosts_[2] + boundsWithoutGhosts_[3]) / 2;
        localGlobalBounds_[4 + i].z = boundsWithoutGhosts_[4 + i];
      }
    }

    // The output extent is computed based on the previous results
    // It is inflated or deflated for boundaries that will contain
    // periodic ghosts
    for(int i = 0; i < 3; i++) {
      if(!(localGlobalBounds_[2 * i].isBound == 1
           && localGlobalBounds_[2 * i + 1].isBound == 1)) {
        outExtent_[2 * i] = inExtent_[2 * i] - 1;
        outExtent_[2 * i + 1] = inExtent_[2 * i + 1] + 1;
      } else {
        outExtent_[2 * i] = inExtent_[2 * i];
        outExtent_[2 * i + 1] = inExtent_[2 * i + 1];
      }
    }

    // We compute here if the global boundary has periodic ghosts. In case a
    // process has both a low and high (e.g. w low and x high) global boundary,
    // the boundary does not have periodic ghosts
    for(int i = 0; i < 3; i++) {
      if(localGlobalBounds_[2 * i].isBound == 1
         && localGlobalBounds_[2 * i].isBound
              == localGlobalBounds_[2 * i + 1].isBound) {
        isBoundaryPeriodic_[2 * i] = 0;
        isBoundaryPeriodic_[2 * i + 1] = 0;
      } else {
        isBoundaryPeriodic_[2 * i] = localGlobalBounds_[2 * i].isBound;
        isBoundaryPeriodic_[2 * i + 1] = localGlobalBounds_[2 * i + 1].isBound;
      }
    }
    isOutputExtentComputed_ = true;
  }
  return 1;
};

template <int matchesSize, int metaDataSize>
int ttkPeriodicGhostsGeneration::MarshalAndSendRecv(
  vtkImageData *imageIn,
  std::vector<std::vector<vtkSmartPointer<vtkCharArray>>> &charArrayBoundaries,
  std::vector<std::vector<std::array<ttk::SimplexId, metaDataSize>>>
    &charArrayBoundariesMetaData,
  std::vector<std::array<ttk::SimplexId, matchesSize>> &matches,
  std::vector<vtkSmartPointer<vtkCharArray>> &charArrayBoundariesReceived,
  std::vector<std::array<ttk::SimplexId, metaDataSize>>
    &charArrayBoundariesMetaDataReceived,
  int dim) {

  int *default_VOI = imageIn->GetExtent();
  std::array<ttk::SimplexId, 6> VOI;
  ttk::SimplexId matchesNumber = matches.size();
  std::vector<std::vector<int>> additionalNeighbors(
    ttk::MPIsize_, std::vector<int>{});
  for(int i = 0; i < matchesNumber; i++) {
    // Extent of the vtkImageData to extract
    VOI = {default_VOI[0], default_VOI[1], default_VOI[2],
           default_VOI[3], default_VOI[4], default_VOI[5]};
    // It is modified based on the direction of the match
    for(int k = 1; k <= dim; k++) {
      if(matches[i][k] % 2 == 0) {
        VOI[matches[i][k] + 1] = VOI[matches[i][k]];
      } else {
        VOI[matches[i][k] - 1] = VOI[matches[i][k]];
      }
    }
    // The intersections between the extracted extent and the neighbors
    // are computed to keep track of the new neighbors
    for(const auto neigh : neighbors_) {
      if(neigh != matches[i][0]) {
        if(!(neighborVertexBBoxes_[neigh][0] == 0
             && neighborVertexBBoxes_[neigh][1] == 0
             && neighborVertexBBoxes_[neigh][2] == 0
             && neighborVertexBBoxes_[neigh][3] == 0
             && neighborVertexBBoxes_[neigh][4] == 0
             && neighborVertexBBoxes_[neigh][5] == 0)) {
          if(VOI[0] - inExtent_[0] <= neighborVertexBBoxes_[neigh][1]
             && VOI[1] - inExtent_[0] >= neighborVertexBBoxes_[neigh][0]
             && VOI[2] - inExtent_[2] <= neighborVertexBBoxes_[neigh][3]
             && VOI[3] - inExtent_[2] >= neighborVertexBBoxes_[neigh][2]
             && VOI[4] - inExtent_[4] <= neighborVertexBBoxes_[neigh][5]
             && VOI[5] - inExtent_[4] >= neighborVertexBBoxes_[neigh][4]) {
            if(std::find(additionalNeighbors[matches[i][0]].begin(),
                         additionalNeighbors[matches[i][0]].end(), neigh)
               == additionalNeighbors[matches[i][0]].end()) {
              additionalNeighbors[matches[i][0]].push_back(neigh);
            }
          }
        }
      }
    }
    // The vtkImageData is extracted
    vtkSmartPointer<vtkExtractVOI> extractVOI
      = vtkSmartPointer<vtkExtractVOI>::New();
    extractVOI->SetInputData(imageIn);
    extractVOI->SetVOI(VOI[0], VOI[1], VOI[2], VOI[3], VOI[4], VOI[5]);
    extractVOI->Update();
    vtkSmartPointer<vtkImageData> extracted = extractVOI->GetOutput();
    vtkSmartPointer<vtkCharArray> buffer = vtkSmartPointer<vtkCharArray>::New();
    // The vtkImageData is converted to a vtkCharArray (for sending)
    if(vtkCommunicator::MarshalDataObject(extracted, buffer) == 0) {
      printErr("Marshalling failed!");
    };
    charArrayBoundaries[matches[i][0]].emplace_back(buffer);
    std::array<ttk::SimplexId, metaDataSize> metaData;
    for(int j = 0; j < metaDataSize; j++) {
      metaData.at(j) = matches[i][j + 1];
    }
    // The type of match is stored as meta data
    charArrayBoundariesMetaData[matches[i][0]].emplace_back(metaData);
  }

  ttk::SimplexId recv_size{0};
  ttk::SimplexId send_size{0};
  std::array<ttk::SimplexId, metaDataSize + 1> sendMetadata;
  std::array<ttk::SimplexId, metaDataSize + 1> recvMetaData;
  // All the extracted image data are exchanged between processes along with
  // their meta data and their new neighbors
  for(int i = 0; i < ttk::MPIsize_; i++) {
    for(int j = 0; j < static_cast<int>(charArrayBoundaries[i].size()); j++) {
      send_size = charArrayBoundaries[i][j]->GetNumberOfTuples();
      MPI_Sendrecv(&send_size, 1, ttk::getMPIType(send_size), i, 0, &recv_size,
                   1, ttk::getMPIType(recv_size), i, 0, ttk::MPIcomm_,
                   MPI_STATUS_IGNORE);
      vtkSmartPointer<vtkCharArray> buffer
        = vtkSmartPointer<vtkCharArray>::New();
      buffer->SetNumberOfTuples(recv_size);
      buffer->SetNumberOfComponents(1);
      char *sendPointer = ttkUtils::GetPointer<char>(charArrayBoundaries[i][j]);
      char *recvPointer = ttkUtils::GetPointer<char>(buffer);
      MPI_Sendrecv(sendPointer, send_size, MPI_CHAR, i, 0, recvPointer,
                   recv_size, MPI_CHAR, i, 0, ttk::MPIcomm_, MPI_STATUS_IGNORE);
      charArrayBoundariesReceived.emplace_back(buffer);
      std::copy(charArrayBoundariesMetaData[i][j].begin(),
                charArrayBoundariesMetaData[i][j].begin() + metaDataSize,
                sendMetadata.begin());
      sendMetadata.at(metaDataSize) = additionalNeighbors[i].size();
      MPI_Sendrecv(sendMetadata.data(), metaDataSize + 1,
                   ttk::getMPIType(recv_size), i, 0, recvMetaData.data(),
                   metaDataSize + 1, ttk::getMPIType(recv_size), i, 0,
                   ttk::MPIcomm_, MPI_STATUS_IGNORE);
      std::array<ttk::SimplexId, metaDataSize> aux;
      std::copy(recvMetaData.begin(), recvMetaData.end() - 1, aux.begin());
      charArrayBoundariesMetaDataReceived.emplace_back(aux);
      std::vector<int> receivedNeighbors(recvMetaData.back());
      MPI_Sendrecv(additionalNeighbors[i].data(), additionalNeighbors[i].size(),
                   MPI_INTEGER, i, 0, receivedNeighbors.data(),
                   recvMetaData.back(), MPI_INTEGER, i, 0, ttk::MPIcomm_,
                   MPI_STATUS_IGNORE);
      for(const int neigh : receivedNeighbors) {
        if(std::find(neighbors_.begin(), neighbors_.end(), neigh)
           == neighbors_.end()) {
          neighbors_.push_back(neigh);
        }
      }
    }
  }
  return 1;
}

template <typename boundaryType>
int ttkPeriodicGhostsGeneration::UnMarshalAndMerge(
  std::vector<boundaryType> &metaDataReceived,
  std::vector<vtkSmartPointer<vtkCharArray>> &boundariesReceived,
  boundaryType direction,
  int mergeDirection,
  vtkImageData *mergedImage) {
  // Search for the right received extracted image data using the meta
  // data and the boundary direction
  auto it
    = std::find(metaDataReceived.begin(), metaDataReceived.end(), direction);
  if(it != metaDataReceived.end()) {
    vtkNew<vtkStructuredPoints> id;
    vtkNew<vtkImageData> aux;
    // If the right extracted image is found, it is converted back to a
    // vtkImageData
    if(vtkCommunicator::UnMarshalDataObject(
         boundariesReceived[std::distance(metaDataReceived.begin(), it)], id)
       == 0) {
      printErr("UnMarshaling failed!");
      return 0;
    };
    // Merged into mergedImage
    this->MergeImageAndSlice(mergedImage, id, aux, mergeDirection);
    mergedImage->DeepCopy(aux);
  }
  return 1;
}

template <typename boundaryType>
int ttkPeriodicGhostsGeneration::UnMarshalAndCopy(
  std::vector<boundaryType> &metaDataReceived,
  std::vector<vtkSmartPointer<vtkCharArray>> &boundariesReceived,
  boundaryType direction,
  vtkImageData *mergedImage) {
  // Search for the right received extracted image data using the meta
  // data and the boundary direction
  auto it
    = std::find(metaDataReceived.begin(), metaDataReceived.end(), direction);
  if(it != metaDataReceived.end()) {
    vtkNew<vtkStructuredPoints> id;
    // If the right extracted image is found, it is converted back to a
    // vtkImageData
    if(vtkCommunicator::UnMarshalDataObject(
         boundariesReceived[std::distance(metaDataReceived.begin(), it)], id)
       == 0) {
      printErr("UnMarshaling failed!");
      return 0;
    };
    // Copied into mergedImage
    mergedImage->DeepCopy(id);
  }
  return 1;
}

int ttkPeriodicGhostsGeneration::MergeDataArrays(
  vtkDataArray *imageArray,
  vtkDataArray *sliceArray,
  vtkSmartPointer<vtkDataArray> &currentArray,
  int direction,
  int dims[3],
  unsigned char ghostValue,
  ttk::SimplexId numberOfSimplices,
  ttk::SimplexId numberOfTuples) {
  std::string arrayName(imageArray->GetName());

#ifndef TTK_ENABLE_KAMIKAZE
  if(!sliceArray) {
    printErr("Array " + arrayName
             + " is not present in the Data of the second vtkImageData");
    return 0;
  }
#endif
  currentArray->SetNumberOfComponents(1);
  currentArray->SetNumberOfTuples(numberOfSimplices);
  currentArray->SetName(arrayName.c_str());
  // Special treatment for ghost arrays
  if(std::strcmp(currentArray->GetName(), "vtkGhostType") == 0) {
    if(ghostValue != 0) {
      sliceArray->SetNumberOfTuples(numberOfTuples);
      sliceArray->Fill(ghostValue);
    } else {
      for(int i = 0; i < numberOfTuples; i++) {
        if(sliceArray->GetTuple1(i) == 1) {
          sliceArray->SetTuple1(i, 16);
        }
      }
    }
  }
  int sliceCounter = 0;
  int imageCounter = 0;
  int counter = 0;
  // We define the function addSlice based on the direction of the merge
  std::function<bool(int, int, int, int[3])> addSlice;
  switch(direction) {
    case 0:
      addSlice = [](int x, int ttkNotUsed(y), int ttkNotUsed(z),
                    int ttkNotUsed(dimensions)[3]) { return x == 0; };
      break;
    case 1:
      addSlice = [](int x, int ttkNotUsed(y), int ttkNotUsed(z),
                    int dimensions[3]) { return x == dimensions[0] - 1; };
      break;
    case 2:
      addSlice = [](int ttkNotUsed(x), int y, int ttkNotUsed(z),
                    int ttkNotUsed(dimensions)[3]) { return y == 0; };
      break;
    case 3:
      addSlice = [](int ttkNotUsed(x), int y, int ttkNotUsed(z),
                    int dimensions[3]) { return y == dimensions[1] - 1; };
      break;
    case 4:
      addSlice = [](int ttkNotUsed(x), int ttkNotUsed(y), int z,
                    int ttkNotUsed(dimensions)[3]) { return z == 0; };
      break;
    case 5:
      addSlice = [](int ttkNotUsed(x), int ttkNotUsed(y), int z,
                    int dimensions[3]) { return z == dimensions[2] - 1; };
      break;
  }

  // Do the actual merging
  for(int z = 0; z < dims[2]; z++) {
    for(int y = 0; y < dims[1]; y++) {
      for(int x = 0; x < dims[0]; x++) {
        if(addSlice(x, y, z, dims)) {
          currentArray->SetTuple1(counter, sliceArray->GetTuple1(sliceCounter));
          sliceCounter++;
        } else {
          currentArray->SetTuple1(counter, imageArray->GetTuple1(imageCounter));
          imageCounter++;
        }
        counter++;
      }
    }
  }
  return 1;
}

int ttkPeriodicGhostsGeneration::MergeImageAndSlice(vtkImageData *image,
                                                    vtkImageData *slice,
                                                    vtkImageData *mergedImage,
                                                    int direction) {
  vtkDataArray *arrayImage = image->GetCellData()->GetArray("Cell Type");
  vtkDataArray *arraySlice = slice->GetCellData()->GetArray("Cell Type");
  if(!(arrayImage && arraySlice) && arrayImage) {
    image->GetCellData()->RemoveArray("Cell Type");
  }
  if(!(arrayImage && arraySlice) && arraySlice) {
    slice->GetCellData()->RemoveArray("Cell Type");
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(image->GetPointData()->GetNumberOfArrays()
     != slice->GetPointData()->GetNumberOfArrays()) {
    printErr("The two vtkImageData objects to merge don't have the same number "
             "of arrays for their Point Data");
    return 0;
  }

  if(image->GetCellData()->GetNumberOfArrays()
     != slice->GetCellData()->GetNumberOfArrays()) {
    printErr("The two vtkImageData objects to merge don't have the same number "
             "of arrays for their Cell Data");
    return 0;
  }
#endif
  mergedImage->SetSpacing(image->GetSpacing());
  mergedImage->Initialize();
  int extentImage[6];
  image->GetExtent(extentImage);
  if(direction % 2 == 0) {
    extentImage[direction] -= 1;
  } else {
    extentImage[direction] += 1;
  }
  mergedImage->SetExtent(extentImage);
  int dims[3];
  mergedImage->GetDimensions(dims);
  int cellDims[3] = {std::max(dims[0] - 1, 1), std::max(dims[1] - 1, 1),
                     std::max(dims[2] - 1, 1)};
  int numberOfPoints = mergedImage->GetNumberOfPoints();
  int imageDims[3];
  image->GetDimensions(imageDims);
  int sliceDims[3];
  slice->GetDimensions(sliceDims);

  for(int array = 0; array < image->GetPointData()->GetNumberOfArrays();
      array++) {
    vtkDataArray *imageArray = image->GetPointData()->GetArray(array);
    vtkDataArray *sliceArray
      = slice->GetPointData()->GetArray(imageArray->GetName());
    vtkSmartPointer<vtkDataArray> currentArray
      = vtkSmartPointer<vtkDataArray>::Take(imageArray->NewInstance());
    this->MergeDataArrays(imageArray, sliceArray, currentArray, direction, dims,
                          vtkDataSetAttributes::DUPLICATEPOINT, numberOfPoints,
                          slice->GetNumberOfPoints());
    mergedImage->GetPointData()->AddArray(currentArray);
  }
  int numberOfCells = mergedImage->GetNumberOfCells();

  for(int array = 0; array < image->GetCellData()->GetNumberOfArrays();
      array++) {
    vtkDataArray *imageArray = image->GetCellData()->GetArray(array);
    std::string arrayName(imageArray->GetName());
    vtkDataArray *sliceArray
      = slice->GetCellData()->GetArray(arrayName.c_str());
    if(std::strcmp(arrayName.c_str(), "Cell Type") == 0) {
      continue;
    }
    vtkSmartPointer<vtkDataArray> currentArray
      = vtkSmartPointer<vtkDataArray>::Take(imageArray->NewInstance());
    // If the direction is low, the new cells created belong to the process
    if(direction == 0 || direction == 2 || direction == 4) {
      this->MergeDataArrays(imageArray, sliceArray, currentArray, direction,
                            cellDims, 0, numberOfCells,
                            slice->GetNumberOfCells());
    } else {
      // The value of ghost is not 1, as it would make the cell disappear in
      // Paraview
      this->MergeDataArrays(imageArray, sliceArray, currentArray, direction,
                            cellDims, vtkDataSetAttributes::EXTERIORCELL,
                            numberOfCells, slice->GetNumberOfCells());
    }

    mergedImage->GetCellData()->AddArray(currentArray);
  }

  vtkNew<vtkCharArray> cellTypes;
  cellTypes->SetName("Cell Type");
  cellTypes->SetNumberOfTuples(numberOfCells);
  cellTypes->Fill(image->GetCellType(0));
  mergedImage->GetCellData()->AddArray(cellTypes);

  return 1;
};

int ttkPeriodicGhostsGeneration::MPIPeriodicGhostPipelinePreconditioning(
  vtkImageData *imageIn, vtkImageData *imageOut) {

  if(!ttk::isRunningWithMPI()) {
    return 0;
  }

  auto other = [](ttk::SimplexId i) {
    if(i % 2 == 1) {
      return i - 1;
    }
    return i + 1;
  };

  int dimensionality = imageIn->GetDataDimension();

  ttk::RegularGridTriangulation *triangulation
    = (ttk::RegularGridTriangulation *)GetTriangulation(imageIn);
  neighborVertexBBoxes_ = triangulation->getNeighborVertexBBoxes();

  // Computation of the order in which the merging of the ghosts should be done
  // The first dimension is the longest
  std::map<int, ttk::SimplexId> dimensionOrder;
  dimensionOrder[0] = inExtent_[1] - inExtent_[0];
  dimensionOrder[2] = inExtent_[3] - inExtent_[2];
  dimensionOrder[4] = inExtent_[5] - inExtent_[4];

  int firstDim{-1};
  int secondDim{-1};
  int thirdDim{-1};

  auto maxMap = [](const std::pair<int, ttk::SimplexId> &p1,
                   const std::pair<int, ttk::SimplexId> &p2) {
    return p1.second <= p2.second;
  };
  auto pr = std::max_element(
    std::begin(dimensionOrder), std::end(dimensionOrder), maxMap);
  firstDim = pr->first;
  dimensionOrder[pr->first] = std::round(pr->second / 2);
  pr = std::max_element(
    std::begin(dimensionOrder), std::end(dimensionOrder), maxMap);
  if(firstDim != pr->first) {
    secondDim = pr->first;
  }
  dimensionOrder[pr->first] = std::round(pr->second / 2);

  pr = std::max_element(
    std::begin(dimensionOrder), std::end(dimensionOrder), maxMap);
  if(secondDim == -1) {
    secondDim = pr->first;
  }
  thirdDim = 6 - (firstDim + secondDim);

  this->ComputeOutputExtent();

  // Preparation of the MPI type for the matches
  MPI_Datatype partialGlobalBoundMPI;
  std::vector<ttk::periodicGhosts::partialGlobalBound> allLocalGlobalBounds(
    ttk::MPIsize_ * 6, ttk::periodicGhosts::partialGlobalBound{});
  MPI_Datatype types[]
    = {MPI_UNSIGNED_CHAR, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  int lengths[] = {1, 1, 1, 1};
  const long int mpi_offsets[]
    = {offsetof(ttk::periodicGhosts::partialGlobalBound, isBound),
       offsetof(ttk::periodicGhosts::partialGlobalBound, x),
       offsetof(ttk::periodicGhosts::partialGlobalBound, y),
       offsetof(ttk::periodicGhosts::partialGlobalBound, z)};
  MPI_Type_create_struct(
    4, lengths, mpi_offsets, types, &partialGlobalBoundMPI);
  MPI_Type_commit(&partialGlobalBoundMPI);

  MPI_Allgather(localGlobalBounds_.data(), 6, partialGlobalBoundMPI,
                allLocalGlobalBounds.data(), 6, partialGlobalBoundMPI,
                ttk::MPIcomm_);

  // Computation of the 1D matches
  std::vector<std::array<ttk::SimplexId, 3>> matches;
  for(int i = 0; i < ttk::MPIsize_; i++) {
    if(i != ttk::MPIrank_) {
      for(int j = 0; j < 6; j++) {
        bool isIn = false;
        if(!(localGlobalBounds_[other(j)].isBound != 0
             && localGlobalBounds_[j].isBound != 0)) {
          if(localGlobalBounds_[other(j)].isBound != 0
             && allLocalGlobalBounds[i * 6 + j].isBound != 0) {
            if(0 <= j && j <= 1) {
              isIn
                = (boundsWithoutGhosts_[2] <= allLocalGlobalBounds[i * 6 + j].y
                   && boundsWithoutGhosts_[3]
                        >= allLocalGlobalBounds[i * 6 + j].y)
                  && (boundsWithoutGhosts_[4]
                        <= allLocalGlobalBounds[i * 6 + j].z
                      && boundsWithoutGhosts_[5]
                           >= allLocalGlobalBounds[i * 6 + j].z);
            }
            if(2 <= j && j <= 3) {
              isIn
                = (boundsWithoutGhosts_[0] <= allLocalGlobalBounds[i * 6 + j].x
                   && boundsWithoutGhosts_[1]
                        >= allLocalGlobalBounds[i * 6 + j].x)
                  && (boundsWithoutGhosts_[4]
                        <= allLocalGlobalBounds[i * 6 + j].z
                      && boundsWithoutGhosts_[5]
                           >= allLocalGlobalBounds[i * 6 + j].z);
            }
            if(4 <= j && j <= 5) {
              isIn
                = (boundsWithoutGhosts_[0] <= allLocalGlobalBounds[i * 6 + j].x
                   && boundsWithoutGhosts_[1]
                        >= allLocalGlobalBounds[i * 6 + j].x)
                  && (boundsWithoutGhosts_[2]
                        <= allLocalGlobalBounds[i * 6 + j].y
                      && boundsWithoutGhosts_[3]
                           >= allLocalGlobalBounds[i * 6 + j].y);
            }
            if(isIn) {
              matches.emplace_back(
                std::array<ttk::SimplexId, 3>{i, other(j), j});
              if(std::find(neighbors_.begin(), neighbors_.end(), i)
                 == neighbors_.end()) {
                neighbors_.push_back(i);
              }
            }
          }
        }
      }
    }
  }
  // Computation of the 2D matches
  std::vector<std::array<ttk::SimplexId, 2>> local2DBounds;
  std::vector<std::array<ttk::SimplexId, 3>> matches_2D;
  if(dimensionality >= 2) {
    for(int i = 0; i < 4; i++) {
      for(int j = i + 1; j < 6; j++) {
        if((abs(i - j) == 1 && i % 2 == 1) || abs(i - j) >= 2) {
          if((localGlobalBounds_[i].isBound != 0
              && localGlobalBounds_[j].isBound != 0)
             && !(localGlobalBounds_[other(i)].isBound != 0
                  && localGlobalBounds_[other(j)].isBound != 0)) {
            local2DBounds.emplace_back(std::array<ttk::SimplexId, 2>{i, j});
          }
        }
      }
    }

    for(int i = 0; i < static_cast<int>(local2DBounds.size()); i++) {
      for(int j = 0; j < ttk::MPIsize_; j++) {
        if(j != ttk::MPIrank_) {
          bool isIn = false;
          if((allLocalGlobalBounds[j * 6 + other(local2DBounds[i][0])].isBound
              != 0)
             && (allLocalGlobalBounds[j * 6 + other(local2DBounds[i][1])]
                   .isBound
                 != 0)) {
            ttk::SimplexId dirs[2]
              = {other(local2DBounds[i][0]), other(local2DBounds[i][1])};
            std::sort(dirs, dirs + 2);
            if((dirs[0] < 2 && dirs[1] >= 2 && dirs[1] < 4)) {
              isIn = (boundsWithoutGhosts_[4]
                        <= allLocalGlobalBounds[j * 6 + dirs[0]].z
                      && boundsWithoutGhosts_[5]
                           >= allLocalGlobalBounds[j * 6 + dirs[0]].z);
            }
            if((dirs[0] < 2 && dirs[1] >= 4)) {
              isIn = (boundsWithoutGhosts_[2]
                        <= allLocalGlobalBounds[j * 6 + dirs[0]].y
                      && boundsWithoutGhosts_[3]
                           >= allLocalGlobalBounds[j * 6 + dirs[0]].y);
            }
            if((dirs[0] >= 2 && dirs[0] < 4 && dirs[1] >= 4)) {
              isIn = (boundsWithoutGhosts_[0]
                        <= allLocalGlobalBounds[j * 6 + dirs[0]].x
                      && boundsWithoutGhosts_[1]
                           >= allLocalGlobalBounds[j * 6 + dirs[0]].x);
            }
            if(isIn) {
              matches_2D.emplace_back(std::array<ttk::SimplexId, 3>{
                j, local2DBounds[i][0], local2DBounds[i][1]});
              if(std::find(neighbors_.begin(), neighbors_.end(), j)
                 == neighbors_.end()) {
                neighbors_.push_back(j);
              }
            }
          }
        }
      }
    }
  }
  // Computation of the 3D matches
  std::vector<std::array<ttk::SimplexId, 3>> local3DBounds;
  std::vector<std::array<ttk::SimplexId, 4>> matches_3D;
  if(dimensionality == 3) {
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++) {
          if(localGlobalBounds_[i].isBound != 0
             && localGlobalBounds_[2 + j].isBound != 0
             && localGlobalBounds_[4 + k].isBound != 0) {
            local3DBounds.emplace_back(
              std::array<ttk::SimplexId, 3>{i, 2 + j, 4 + k});
          }
        }
      }
    }
    for(int i = 0; i < static_cast<int>(local3DBounds.size()); i++) {
      for(int j = 0; j < ttk::MPIsize_; j++) {
        if(j != ttk::MPIrank_) {
          if((allLocalGlobalBounds[j * 6 + other(local3DBounds[i][0])].isBound
              != 0)
             && (allLocalGlobalBounds[j * 6 + other(local3DBounds[i][1])]
                   .isBound
                 != 0)
             && (allLocalGlobalBounds[j * 6 + other(local3DBounds[i][2])]
                   .isBound
                 != 0)) {
            matches_3D.emplace_back(std::array<ttk::SimplexId, 4>{
              j, local3DBounds[i][0], local3DBounds[i][1],
              local3DBounds[i][2]});
            if(std::find(neighbors_.begin(), neighbors_.end(), j)
               == neighbors_.end()) {
              neighbors_.push_back(j);
            }
          }
        }
      }
    }
  }
  // Now, extract ImageData for 1D boundaries
  std::vector<std::vector<vtkSmartPointer<vtkCharArray>>> charArrayBoundaries(
    ttk::MPIsize_);
  std::vector<std::vector<std::array<ttk::SimplexId, 1>>>
    charArrayBoundariesMetaData(ttk::MPIsize_);
  std::vector<vtkSmartPointer<vtkCharArray>> charArrayBoundariesReceived;
  std::vector<std::array<ttk::SimplexId, 1>>
    charArrayBoundariesMetaDataReceived;
  if(this->MarshalAndSendRecv<3, 1>(
       imageIn, charArrayBoundaries, charArrayBoundariesMetaData, matches,
       charArrayBoundariesReceived, charArrayBoundariesMetaDataReceived, 1)
     == 0) {
    printErr("Error occurred when marshalling messages");
    return 0;
  }

  // Extract 2D boundaries
  std::vector<std::vector<vtkSmartPointer<vtkCharArray>>> charArray2DBoundaries(
    ttk::MPIsize_);
  std::vector<std::vector<std::array<ttk::SimplexId, 2>>>
    charArray2DBoundariesMetaData(ttk::MPIsize_);
  std::vector<vtkSmartPointer<vtkCharArray>> charArray2DBoundariesReceived;
  std::vector<std::array<ttk::SimplexId, 2>>
    charArray2DBoundariesMetaDataReceived;
  if(dimensionality >= 2) {
    if(this->MarshalAndSendRecv<3, 2>(imageIn, charArray2DBoundaries,
                                      charArray2DBoundariesMetaData, matches_2D,
                                      charArray2DBoundariesReceived,
                                      charArray2DBoundariesMetaDataReceived, 2)
       == 0) {
      return 0;
    }
  }
  // Now, same for 3D boundaries
  std::vector<std::vector<vtkSmartPointer<vtkCharArray>>> charArray3DBoundaries(
    ttk::MPIsize_);
  std::vector<std::vector<std::array<ttk::SimplexId, 3>>>
    charArray3DBoundariesMetaData(ttk::MPIsize_);
  std::vector<vtkSmartPointer<vtkCharArray>> charArray3DBoundariesReceived;
  std::vector<std::array<ttk::SimplexId, 3>>
    charArray3DBoundariesMetaDataReceived;
  if(dimensionality == 3) {
    if(this->MarshalAndSendRecv<4, 3>(imageIn, charArray3DBoundaries,
                                      charArray3DBoundariesMetaData, matches_3D,
                                      charArray3DBoundariesReceived,
                                      charArray3DBoundariesMetaDataReceived, 3)
       == 0) {
      return 0;
    }
  }
  imageOut->DeepCopy(imageIn);

  // Merge in the first direction (low and high)
  for(int dir = firstDim; dir < firstDim + 2; dir++) {
    if(this->UnMarshalAndMerge<std::array<ttk::SimplexId, 1>>(
         charArrayBoundariesMetaDataReceived, charArrayBoundariesReceived,
         std::array<ttk::SimplexId, 1>{other(dir)}, dir, imageOut)
       == 0) {
      return 0;
    }
  }
  if(dimensionality >= 2) {
    // Merge in the second direction
    for(int dir = secondDim; dir < secondDim + 2; dir++) {
      vtkSmartPointer<vtkImageData> mergedImage
        = vtkSmartPointer<vtkImageData>::New();
      if(this->UnMarshalAndCopy<std::array<ttk::SimplexId, 1>>(
           charArrayBoundariesMetaDataReceived, charArrayBoundariesReceived,
           std::array<ttk::SimplexId, 1>{other(dir)}, mergedImage)
         == 0) {
        return 0;
      }
      if(mergedImage->GetNumberOfPoints() > 0) {
        for(int dir_2D = firstDim; dir_2D < firstDim + 2; dir_2D++) {
          std::array<ttk::SimplexId, 2> merge_dir{
            std::min(other(dir), other(dir_2D)),
            std::max(other(dir_2D), other(dir))};
          if(this->UnMarshalAndMerge<std::array<ttk::SimplexId, 2>>(
               charArray2DBoundariesMetaDataReceived,
               charArray2DBoundariesReceived, merge_dir, dir_2D, mergedImage)
             == 0) {
            return 0;
          }
        }
        vtkSmartPointer<vtkImageData> aux
          = vtkSmartPointer<vtkImageData>::New();
        if(this->MergeImageAndSlice(imageOut, mergedImage, aux, dir) == 0) {
          return 0;
        }
        imageOut->DeepCopy(aux);
      }
    }
  }
  if(dimensionality == 3) {
    // Merge in the third direction
    for(int dir = thirdDim; dir < thirdDim + 2; dir++) {
      vtkNew<vtkImageData> mergedImage1;
      if(this->UnMarshalAndCopy<std::array<ttk::SimplexId, 1>>(
           charArrayBoundariesMetaDataReceived, charArrayBoundariesReceived,
           std::array<ttk::SimplexId, 1>{other(dir)}, mergedImage1)
         == 0) {
        return 0;
      }
      if(mergedImage1->GetNumberOfPoints() > 0) {
        for(int dir_2D = firstDim; dir_2D < firstDim + 2; dir_2D++) {
          std::array<ttk::SimplexId, 2> merge_dir{
            std::min(other(dir), other(dir_2D)),
            std::max(other(dir_2D), other(dir))};
          if(this->UnMarshalAndMerge<std::array<ttk::SimplexId, 2>>(
               charArray2DBoundariesMetaDataReceived,
               charArray2DBoundariesReceived, merge_dir, dir_2D, mergedImage1)
             == 0) {
            return 0;
          }
        }
        for(int dir_2D = secondDim; dir_2D < secondDim + 2; dir_2D++) {
          vtkNew<vtkImageData> mergedImage2;
          std::array<ttk::SimplexId, 2> merge_dir_2D{
            std::min(other(dir), other(dir_2D)),
            std::max(other(dir_2D), other(dir))};
          if(this->UnMarshalAndCopy<std::array<ttk::SimplexId, 2>>(
               charArray2DBoundariesMetaDataReceived,
               charArray2DBoundariesReceived, merge_dir_2D, mergedImage2)
             == 0) {
            return 0;
          }
          if(mergedImage2->GetNumberOfPoints() > 0) {
            for(int dir_3D = firstDim; dir_3D < firstDim + 2; dir_3D++) {
              std::array<ttk::SimplexId, 3> merge_dir_3D{
                other(dir), other(dir_2D), other(dir_3D)};
              std::sort(merge_dir_3D.begin(), merge_dir_3D.end());
              if(this->UnMarshalAndMerge<std::array<ttk::SimplexId, 3>>(
                   charArray3DBoundariesMetaDataReceived,
                   charArray3DBoundariesReceived, merge_dir_3D, dir_3D,
                   mergedImage2)
                 == 0) {
                return 0;
              }
            }
            vtkNew<vtkImageData> aux;
            if(this->MergeImageAndSlice(mergedImage1, mergedImage2, aux, dir_2D)
               == 0) {
              return 0;
            }
            mergedImage1->DeepCopy(aux);
          }
        }
        if(mergedImage1->GetNumberOfPoints() > 0) {
          vtkNew<vtkImageData> aux;
          if(this->MergeImageAndSlice(imageOut, mergedImage1, aux, dir) == 0) {
            return 0;
          }
          imageOut->DeepCopy(aux);
        }
      }
    }
  }

  int nbArrays = imageOut->GetPointData()->GetNumberOfArrays();
  for(int i = nbArrays - 1; i >= 0; i--) {
    vtkDataArray *array = imageOut->GetPointData()->GetArray(i);
    std::string arrayName = std::string(array->GetName());
    if(arrayName.rfind("_Order") == (arrayName.size() - 6)) {
      imageOut->GetPointData()->RemoveArray(i);
    }
  }
  return 1;
};

#endif
int ttkPeriodicGhostsGeneration::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkImageData *imageIn = vtkImageData::GetData(inputVector[0]);
  vtkImageData *imageOut = vtkImageData::GetData(outputVector);
#ifdef TTK_ENABLE_MPI
  this->MPIPeriodicGhostPipelinePreconditioning(imageIn, imageOut);
#else
  imageOut->ShallowCopy(imageIn);
#endif

  // return success
  return 1;
};
