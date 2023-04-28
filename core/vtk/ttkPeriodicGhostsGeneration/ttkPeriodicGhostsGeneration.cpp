#include "vtkStreamingDemandDrivenPipeline.h"
#include <ttkPeriodicGhostsGeneration.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkImageAppend.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkStructuredPoints.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPeriodicGhostsGeneration);

ttkPeriodicGhostsGeneration::ttkPeriodicGhostsGeneration() {
  this->setDebugMsgPrefix("PeriodicGhostsGeneration");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkPeriodicGhostsGeneration::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkPeriodicGhostsGeneration::RequestUpdateExtent(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  vtkImageData *image = vtkImageData::GetData(inputVector[0]);
  this->ComputeOutputExtent(image);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
              outExtent_[0] + 1, outExtent_[1] - 1, outExtent_[2] + 1,
              outExtent_[3] - 1, outExtent_[4] + 1, outExtent_[5] - 1);
  return 1;
}

int ttkPeriodicGhostsGeneration::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVectors,
  vtkInformationVector *outputVector) {
  vtkImageData *image = vtkImageData::GetData(inputVectors[0]);
  this->ComputeOutputExtent(image);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outExtent_.data(), 6);
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
  return 1;
}

int ttkPeriodicGhostsGeneration::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkPeriodicGhostsGeneration::ComputeOutputExtent(vtkDataSet *input) {
  // TODO: fix problem of direction that don't match because on same process
  // TODO: check for 4 processes and 2 processes
  if(!isOutputExtentComputed_) {
    vtkImageData *imageIn;
    if(input->IsA("vtkImageData")) {
      imageIn = vtkImageData::SafeDownCast(input);
    } else {
      printErr("Invalid data input type for periodicTriangulation computation");
      return -1;
    }

    std::array<double, 6> tempGlobalBounds{};
    double bounds[6];
    imageIn->GetBounds(bounds);
    // Reorganize bounds to only execute Allreduce twice
    std::array<double, 6> tempBounds = {
      bounds[0], bounds[2], bounds[4], bounds[1], bounds[3], bounds[5],
    };

    // Compute and send to all processes the lower bounds of the data set
    MPI_Allreduce(tempBounds.data(), tempGlobalBounds.data(), 3, MPI_DOUBLE,
                  MPI_MIN, ttk::MPIcomm_);
    // Compute and send to all processes the higher bounds of the data set
    MPI_Allreduce(&tempBounds[3], &tempGlobalBounds[3], 3, MPI_DOUBLE, MPI_MAX,
                  ttk::MPIcomm_);

    // re-order tempGlobalBounds
    globalBounds_
      = {tempGlobalBounds[0], tempGlobalBounds[3], tempGlobalBounds[1],
         tempGlobalBounds[4], tempGlobalBounds[2], tempGlobalBounds[5]};

    double spacing[3];
    imageIn->GetSpacing(spacing);
    outExtent_[0] = static_cast<int>(round(globalBounds_[0] / spacing[0])) - 1;
    outExtent_[1] = static_cast<int>(round(globalBounds_[1] / spacing[0])) + 1;
    outExtent_[2] = static_cast<int>(round(globalBounds_[2] / spacing[1])) - 1;
    outExtent_[3] = static_cast<int>(round(globalBounds_[3] / spacing[1])) + 1;
    outExtent_[4] = static_cast<int>(round(globalBounds_[4] / spacing[2])) - 1;
    outExtent_[5] = static_cast<int>(round(globalBounds_[5] / spacing[2])) + 1;
    boundsWithoutGhosts_
      = {bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]};
    for(int i = 0; i < 3; i++) {
      if(bounds[2 * i] != globalBounds_[2 * i]) {
        boundsWithoutGhosts_[2 * i] = bounds[2 * i] + spacing[0];
      } else {
      }
      if(bounds[2 * i + 1] != globalBounds_[2 * i + 1]) {
        boundsWithoutGhosts_[2 * i + 1] = bounds[2 * i + 1] - spacing[0];
      } else {
      }
    }
    isOutputExtentComputed_ = true;
  }
  return 1;
};

int ttkPeriodicGhostsGeneration::MergeImageAppendAndSlice(
  vtkImageData *image,
  vtkImageData *slice,
  vtkImageData *mergedImage,
  int direction) {

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
  int numberOfPoints = mergedImage->GetNumberOfPoints();
  for(int array = 0; array < image->GetPointData()->GetNumberOfArrays();
      array++) {
    vtkDataArray *imageArray = image->GetPointData()->GetArray(array);
    std::string arrayName(imageArray->GetName());
    vtkDataArray *sliceArray
      = slice->GetPointData()->GetArray(arrayName.c_str());
#ifndef TTK_ENABLE_KAMIKAZE
    if(!sliceArray) {
      printErr(
        "Array " + arrayName
        + " is not present in the Point Data of the second vtkImageData");
      return 0;
    }
#endif
    vtkSmartPointer<vtkDataArray> currentArray
      = vtkSmartPointer<vtkDataArray>::Take(imageArray->NewInstance());
    currentArray->SetNumberOfComponents(1);
    currentArray->SetNumberOfTuples(numberOfPoints);
    currentArray->SetName(arrayName.c_str());
    if(std::strcmp(currentArray->GetName(), "vtkGhostType") == 0) {
      sliceArray->SetNumberOfTuples(slice->GetNumberOfPoints());
      sliceArray->Fill(vtkDataSetAttributes::DUPLICATEPOINT);
    }
    int sliceCounter = 0;
    int imageCounter = 0;
    int counter = 0;
    switch(direction) {
      case 0:
        for(int z = 0; z < dims[2]; z++) {
          for(int y = 0; y < dims[1]; y++) {
            for(int x = 0; x < dims[0]; x++) {
              if(x == 0) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 1:
        for(int z = 0; z < dims[2]; z++) {
          for(int y = 0; y < dims[1]; y++) {
            for(int x = 0; x < dims[0]; x++) {
              if(x == dims[0] - 1) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 2:
        for(int z = 0; z < dims[2]; z++) {
          for(int y = 0; y < dims[1]; y++) {
            for(int x = 0; x < dims[0]; x++) {
              if(y == 0) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 3:
        for(int z = 0; z < dims[2]; z++) {
          for(int y = 0; y < dims[1]; y++) {
            for(int x = 0; x < dims[0]; x++) {
              if(y == dims[1] - 1) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 4:
        for(int z = 0; z < dims[2]; z++) {
          for(int y = 0; y < dims[1]; y++) {
            for(int x = 0; x < dims[0]; x++) {
              if(z == 0) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 5:
        for(int z = 0; z < dims[2]; z++) {
          for(int y = 0; y < dims[1]; y++) {
            for(int x = 0; x < dims[0]; x++) {
              if(z == dims[2] - 1) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
    };
    mergedImage->GetPointData()->AddArray(currentArray);
  }
  int numberOfCells = mergedImage->GetNumberOfCells();

  for(int array = 0; array < image->GetCellData()->GetNumberOfArrays();
      array++) {
    vtkDataArray *imageArray = image->GetCellData()->GetArray(array);
    std::string arrayName(imageArray->GetName());
    vtkDataArray *sliceArray
      = slice->GetCellData()->GetArray(arrayName.c_str());
#ifndef TTK_ENABLE_KAMIKAZE
    if(!sliceArray) {
      printErr(
        "Array " + arrayName
        + " is not present in the Point Data of the second vtkImageData");
      return 0;
    }
#endif
    vtkSmartPointer<vtkDataArray> currentArray
      = vtkSmartPointer<vtkDataArray>::Take(imageArray->NewInstance());
    currentArray->SetNumberOfComponents(1);
    currentArray->SetNumberOfTuples(numberOfCells);
    currentArray->SetName(arrayName.c_str());
    if(std::strcmp(currentArray->GetName(), "vtkGhostType") == 0) {
      sliceArray->SetNumberOfTuples(slice->GetNumberOfCells());
      sliceArray->Fill(vtkDataSetAttributes::EXTERIORCELL);
    }
    vtkIdType sliceCounter = 0;
    vtkIdType imageCounter = 0;
    vtkIdType counter = 0;
    switch(direction) {
      case 0:
        for(int z = 0; z < dims[2] - 1; z++) {
          for(int y = 0; y < dims[1] - 1; y++) {
            for(int x = 0; x < dims[0] - 1; x++) {
              if(x == 0) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 1:
        for(int z = 0; z < dims[2] - 1; z++) {
          for(int y = 0; y < dims[1] - 1; y++) {
            for(int x = 0; x < dims[0] - 1; x++) {
              if(x == dims[0] - 2) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 2:
        for(int z = 0; z < dims[2] - 1; z++) {
          for(int y = 0; y < dims[1] - 1; y++) {
            for(int x = 0; x < dims[0] - 1; x++) {
              if(y == 0) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 3:
        for(int z = 0; z < dims[2] - 1; z++) {
          for(int y = 0; y < dims[1] - 1; y++) {
            for(int x = 0; x < dims[0] - 1; x++) {
              if(y == dims[1] - 2) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 4:
        for(int z = 0; z < dims[2] - 1; z++) {
          for(int y = 0; y < dims[1] - 1; y++) {
            for(int x = 0; x < dims[0] - 1; x++) {
              if(z == 0) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
      case 5:
        for(int z = 0; z < dims[2] - 1; z++) {
          for(int y = 0; y < dims[1] - 1; y++) {
            for(int x = 0; x < dims[0] - 1; x++) {
              if(z == dims[2] - 2) {
                currentArray->SetTuple1(
                  counter, sliceArray->GetTuple1(sliceCounter));
                sliceCounter++;
              } else {
                currentArray->SetTuple1(
                  counter, imageArray->GetTuple1(imageCounter));
                imageCounter++;
              }
              counter++;
            }
          }
        }
        break;
    };
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
  vtkInformationVector **inputVectors, vtkInformationVector *outputVector) {

  struct partialGlobalBound {
    unsigned char isBound{0};
    double x{0};
    double y{0};
    double z{0};
  };

  auto other = [](ttk::SimplexId i) {
    if(i % 2 == 1) {
      return i - 1;
    }
    return i + 1;
  };

  vtkImageData *imageIn = vtkImageData::GetData(inputVectors[0]);
  vtkImageData *imageOut = vtkImageData::GetData(outputVector);

  imageOut->ShallowCopy(imageIn);

  std::array<partialGlobalBound, 6> localGlobalBounds;
  this->ComputeOutputExtent(imageIn);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outExtent_.data(), 6);
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());

  /*std::string s = "boundsWithoutGhosts: ";
  for(int i = 0; i < 6; i++) {
    s += std::to_string(boundsWithoutGhosts_[i]) + ", ";
  }
  printMsg(s);*/

  for(int i = 0; i < 2; i++) {
    if(globalBounds_[i] == boundsWithoutGhosts_[i]) {
      localGlobalBounds[i].isBound = 1;
      localGlobalBounds[i].x = boundsWithoutGhosts_[i];
      localGlobalBounds[i].y
        = (boundsWithoutGhosts_[2] + boundsWithoutGhosts_[3]) / 2;
      localGlobalBounds[i].z
        = (boundsWithoutGhosts_[4] + boundsWithoutGhosts_[5]) / 2;
    }
  }

  for(int i = 0; i < 2; i++) {
    if(globalBounds_[2 + i] == boundsWithoutGhosts_[2 + i]) {
      localGlobalBounds[2 + i].isBound = 1;
      localGlobalBounds[2 + i].x
        = (boundsWithoutGhosts_[0] + boundsWithoutGhosts_[1]) / 2;
      localGlobalBounds[2 + i].y = boundsWithoutGhosts_[2 + i];
      localGlobalBounds[2 + i].z
        = (boundsWithoutGhosts_[4] + boundsWithoutGhosts_[5]) / 2;
    }
  }

  for(int i = 0; i < 2; i++) {
    if(globalBounds_[4 + i] == boundsWithoutGhosts_[4 + i]) {
      localGlobalBounds[4 + i].isBound = 1;
      localGlobalBounds[4 + i].x
        = (boundsWithoutGhosts_[0] + boundsWithoutGhosts_[1]) / 2;
      localGlobalBounds[4 + i].y
        = (boundsWithoutGhosts_[2] + boundsWithoutGhosts_[3]) / 2;
      localGlobalBounds[4 + i].z = boundsWithoutGhosts_[4 + i];
    }
  }
  /*std::string s1 = "";
  printMsg("Local: ");
  for(int j = 0; j < 6; j++) {
    s1 += "isBound: " + std::to_string(localGlobalBounds[j].isBound)
          + "mid: " + std::to_string(localGlobalBounds[j].x) + ", "
          + std::to_string(localGlobalBounds[j].y) + ", "
          + std::to_string(localGlobalBounds[j].z) + "\n";
  }
  printMsg(s1);*/

  MPI_Datatype partialGlobalBoundMPI;
  std::vector<partialGlobalBound> allLocalGlobalBounds(
    ttk::MPIsize_ * 6, partialGlobalBound{});
  MPI_Datatype types[]
    = {MPI_UNSIGNED_CHAR, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  int lengths[] = {1, 1, 1, 1};
  const long int mpi_offsets[]
    = {offsetof(partialGlobalBound, isBound), offsetof(partialGlobalBound, x),
       offsetof(partialGlobalBound, y), offsetof(partialGlobalBound, z)};
  MPI_Type_create_struct(
    4, lengths, mpi_offsets, types, &partialGlobalBoundMPI);
  MPI_Type_commit(&partialGlobalBoundMPI);

  MPI_Allgather(localGlobalBounds.data(), 6, partialGlobalBoundMPI,
                allLocalGlobalBounds.data(), 6, partialGlobalBoundMPI,
                ttk::MPIcomm_);
  /*printMsg("allLocalGlobalBounds: ");
   for (int i = 0; i<ttk::MPIsize_; i++){
     std::string s = "";
     for (int j = 0; j < 6; j++){
         s += "isBound: "+std::to_string(allLocalGlobalBounds[6*i+j].isBound)
         + "mid: "+std::to_string(allLocalGlobalBounds[6*i+j].x)+", "
         +std::to_string(allLocalGlobalBounds[6*i+j].y)+", "
         +std::to_string(allLocalGlobalBounds[6*i+j].z)+"\n";
     }
     printMsg(s);
   } */

  std::vector<std::array<ttk::SimplexId, 3>> matches;
  for(int i = 0; i < ttk::MPIsize_; i++) {
    if(i != ttk::MPIrank_) {
      for(int j = 0; j < 6; j++) {
        bool isIn = false;
        if(!(localGlobalBounds[other(j)].isBound != 0
             && localGlobalBounds[j].isBound != 0)) {
          if(localGlobalBounds[other(j)].isBound != 0
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
            }
          }
        }
      }
    }
  }

  std::vector<std::array<ttk::SimplexId, 2>> local2DBounds;
  for(int i = 0; i < 4; i++) {
    for(int j = i + 1; j < 6; j++) {
      if((abs(i - j) == 1 && i % 2 == 1) || abs(i - j) >= 2) {
        if((localGlobalBounds[i].isBound != 0
            && localGlobalBounds[j].isBound != 0)
           && !(localGlobalBounds[other(i)].isBound != 0
                && localGlobalBounds[other(j)].isBound != 0)) {
          local2DBounds.emplace_back(std::array<ttk::SimplexId, 2>{i, j});
        }
      }
    }
  }
  /*std::string s = "local2Dbounds :";
  for (int i = 0; i<local2DBounds.size();i++){
    s +=
  "("+std::to_string(local2DBounds[i][0])+","+std::to_string(local2DBounds[i][1])+"),
  ";
  }
  printMsg(s);*/
  std::vector<std::array<ttk::SimplexId, 3>> matches_2D;
  for(int i = 0; i < local2DBounds.size(); i++) {
    for(int j = 0; j < ttk::MPIsize_; j++) {
      if(j != ttk::MPIrank_) {
        bool isIn = false;
        if((allLocalGlobalBounds[j * 6 + other(local2DBounds[i][0])].isBound
            != 0)
           && (allLocalGlobalBounds[j * 6 + other(local2DBounds[i][1])].isBound
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
          }
        }
      }
    }
  }

  std::vector<std::array<ttk::SimplexId, 3>> local3DBounds;
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        if(localGlobalBounds[i].isBound != 0
           && localGlobalBounds[2 + j].isBound != 0
           && localGlobalBounds[4 + k].isBound != 0) {
          local3DBounds.emplace_back(
            std::array<ttk::SimplexId, 3>{i, 2 + j, 4 + k});
        }
      }
    }
  }
  std::vector<std::array<ttk::SimplexId, 4>> matches_3D;
  for(int i = 0; i < local3DBounds.size(); i++) {
    for(int j = 0; j < ttk::MPIsize_; j++) {
      if(j != ttk::MPIrank_) {
        if((allLocalGlobalBounds[j * 6 + other(local3DBounds[i][0])].isBound
            != 0)
           && (allLocalGlobalBounds[j * 6 + other(local3DBounds[i][1])].isBound
               != 0)
           && (allLocalGlobalBounds[j * 6 + other(local3DBounds[i][2])].isBound
               != 0)) {
          matches_3D.emplace_back(std::array<ttk::SimplexId, 4>{
            j, local3DBounds[i][0], local3DBounds[i][1], local3DBounds[i][2]});
        }
      }
    }
  }
  /*printMsg("2D matches: "+std::to_string(matches_2D.size()));
  for (int i = 0; i<matches_2D.size(); i++){
    printMsg("proc,i,j: "+std::to_string(matches_2D[i][0])+",
  "+std::to_string(matches_2D[i][1])+", "+std::to_string(matches_2D[i][2]));
  }
  printMsg("3D matches: " + std::to_string(matches_3D.size()));
  for(int i = 0; i < matches_3D.size(); i++) {
    printMsg("proc,i,j,k: " + std::to_string(matches_3D[i][0]) + ", "
             + std::to_string(matches_3D[i][1]) + ", "
             + std::to_string(matches_3D[i][2]) + ", "
             + std::to_string(matches_3D[i][3]));
  }*/

  // Now, extract ImageData
  int *default_VOI = imageIn->GetExtent();
  std::array<ttk::SimplexId, 6> VOI;
  std::vector<std::vector<vtkSmartPointer<vtkCharArray>>> charArrayBoundaries(
    ttk::MPIsize_);
  std::vector<std::vector<ttk::SimplexId>> charArrayBoundariesMetaData(
    ttk::MPIsize_);
  for(int i = 0; i < matches.size(); i++) {
    vtkNew<vtkImageAppend> imageAppend;
    // imageAppend->PreserveExtentsOn();
    VOI = {default_VOI[0], default_VOI[1], default_VOI[2],
           default_VOI[3], default_VOI[4], default_VOI[5]};
    imageAppend->SetAppendAxis(static_cast<int>(matches[i][1] / 2));
    if(matches[i][1] % 2 == 0) {
      VOI[matches[i][1] + 1] = VOI[matches[i][1]];
    } else {
      VOI[matches[i][1] - 1] = VOI[matches[i][1]];
    }
    vtkSmartPointer<vtkExtractVOI> extractVOI
      = vtkSmartPointer<vtkExtractVOI>::New();
    extractVOI->SetInputData(imageIn);
    extractVOI->SetVOI(VOI[0], VOI[1], VOI[2], VOI[3], VOI[4], VOI[5]);
    extractVOI->Update();
    vtkSmartPointer<vtkImageData> extracted = extractVOI->GetOutput();
    vtkSmartPointer<vtkCharArray> buffer = vtkSmartPointer<vtkCharArray>::New();
    int *extent = extracted->GetExtent();
    /*printMsg("Extent: "+std::to_string(extent[0])+", "
    +std::to_string(extent[1])+", "
    +std::to_string(extent[2])+", "
    +std::to_string(extent[3])+", "
    +std::to_string(extent[4])+", "
    +std::to_string(extent[5]));*/

    // printMsg("Size of extracted ImageData:
    // "+std::to_string(extracted->GetNumberOfPoints()));
    // charArrayBoundaries[matches[i][0]].emplace_back(vtkNew<vtkCharArray>());
    if(vtkCommunicator::MarshalDataObject(extracted, buffer) == 0) {
      printErr("Marshalling failed!");
    };
    charArrayBoundaries[matches[i][0]].emplace_back(buffer);
    charArrayBoundariesMetaData[matches[i][0]].emplace_back(matches[i][1]);
  }
  /*for (int i = 0; i < ttk::MPIsize_; i++){
    for (int j = 0; j < charArrayBoundaries[i].size(); j++){
      printMsg("Size of charArray:
  "+std::to_string(charArrayBoundaries[i][j]->GetNumberOfTuples()));
    }
  }*/
  ttk::SimplexId recv_size;
  ttk::SimplexId send_size;
  ttk::SimplexId sendMetadata;
  ttk::SimplexId recvMetaData;
  std::vector<vtkSmartPointer<vtkCharArray>> charArrayBoundariesReceived;
  std::vector<ttk::SimplexId> charArrayBoundariesMetaDataReceived;
  for(int i = 0; i < ttk::MPIsize_; i++) {
    for(int j = 0; j < charArrayBoundaries[i].size(); j++) {
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
      sendMetadata = charArrayBoundariesMetaData[i][j];
      MPI_Sendrecv(&sendMetadata, 1, ttk::getMPIType(sendMetadata), i, 0,
                   &recvMetaData, 1, ttk::getMPIType(recvMetaData), i, 0,
                   ttk::MPIcomm_, MPI_STATUS_IGNORE);
      charArrayBoundariesMetaDataReceived.emplace_back(recvMetaData);
    }
  }
  /*for (int i = 0; i < charArrayBoundariesReceived.size(); i++){
    printMsg("Size of received charArray:
  "+std::to_string(charArrayBoundariesReceived[i]->GetNumberOfTuples()));
  }*/

  // Now, extract ImageData
  std::vector<std::vector<vtkSmartPointer<vtkCharArray>>> charArray2DBoundaries(
    ttk::MPIsize_);
  std::vector<std::vector<std::array<ttk::SimplexId, 2>>>
    charArray2DBoundariesMetaData(ttk::MPIsize_);
  for(int i = 0; i < matches_2D.size(); i++) {
    VOI = {default_VOI[0], default_VOI[1], default_VOI[2],
           default_VOI[3], default_VOI[4], default_VOI[5]};
    for(int k = 1; k < 3; k++) {
      if(matches_2D[i][k] % 2 == 0) {
        VOI[matches_2D[i][k] + 1] = VOI[matches_2D[i][k]];
      } else {
        VOI[matches_2D[i][k] - 1] = VOI[matches_2D[i][k]];
      }
    }
    vtkNew<vtkExtractVOI> extractVOI;
    extractVOI->SetInputData(imageIn);
    extractVOI->SetVOI(VOI[0], VOI[1], VOI[2], VOI[3], VOI[4], VOI[5]);
    extractVOI->Update();
    vtkImageData *extracted = extractVOI->GetOutput();
    vtkSmartPointer<vtkCharArray> buffer = vtkSmartPointer<vtkCharArray>::New();
    if(vtkCommunicator::MarshalDataObject(extracted, buffer) == 0) {
      printErr("Marshalling failed!");
    };
    charArray2DBoundaries[matches_2D[i][0]].emplace_back(buffer);
    charArray2DBoundariesMetaData[matches_2D[i][0]].emplace_back(
      std::array<ttk::SimplexId, 2>{matches_2D[i][1], matches_2D[i][2]});
  }
  std::array<ttk::SimplexId, 2> sendMetadataArray;
  std::array<ttk::SimplexId, 2> recvMetaDataArray;
  std::vector<vtkSmartPointer<vtkCharArray>> charArray2DBoundariesReceived;
  std::vector<std::array<ttk::SimplexId, 2>>
    charArray2DBoundariesMetaDataReceived;
  for(int i = 0; i < ttk::MPIsize_; i++) {
    for(int j = 0; j < charArray2DBoundaries[i].size(); j++) {
      send_size = charArray2DBoundaries[i][j]->GetNumberOfTuples();
      MPI_Sendrecv(&send_size, 1, ttk::getMPIType(send_size), i, 0, &recv_size,
                   1, ttk::getMPIType(recv_size), i, 0, ttk::MPIcomm_,
                   MPI_STATUS_IGNORE);
      vtkSmartPointer<vtkCharArray> buffer
        = vtkSmartPointer<vtkCharArray>::New();
      buffer->SetNumberOfTuples(recv_size);
      buffer->SetNumberOfComponents(1);
      char *sendPointer
        = ttkUtils::GetPointer<char>(charArray2DBoundaries[i][j]);
      char *recvPointer = ttkUtils::GetPointer<char>(buffer);
      MPI_Sendrecv(sendPointer, send_size, MPI_CHAR, i, 0, recvPointer,
                   recv_size, MPI_CHAR, i, 0, ttk::MPIcomm_, MPI_STATUS_IGNORE);
      charArray2DBoundariesReceived.emplace_back(buffer);
      sendMetadataArray = charArray2DBoundariesMetaData[i][j];
      MPI_Sendrecv(&sendMetadataArray, 2, ttk::getMPIType(sendMetadata), i, 0,
                   &recvMetaDataArray, 2, ttk::getMPIType(recvMetaData), i, 0,
                   ttk::MPIcomm_, MPI_STATUS_IGNORE);
      charArray2DBoundariesMetaDataReceived.emplace_back(recvMetaDataArray);
    }
  }

  /*for (int i = 0; i < charArray2DBoundariesReceived.size(); i++){
    printMsg("Size of received charArray:
  "+std::to_string(charArray2DBoundariesReceived[i]->GetNumberOfTuples()));
  }*/

  // Now, same for 3D boundaries

  std::vector<std::vector<vtkSmartPointer<vtkCharArray>>> charArray3DBoundaries(
    ttk::MPIsize_);
  std::vector<std::vector<std::array<ttk::SimplexId, 3>>>
    charArray3DBoundariesMetaData(ttk::MPIsize_);
  for(int i = 0; i < matches_3D.size(); i++) {
    VOI = {default_VOI[0], default_VOI[1], default_VOI[2],
           default_VOI[3], default_VOI[4], default_VOI[5]};
    for(int k = 1; k < 4; k++) {
      if(matches_3D[i][k] % 2 == 0) {
        VOI[matches_3D[i][k] + 1] = VOI[matches_3D[i][k]] + 1;
      } else {
        VOI[matches_3D[i][k] - 1] = VOI[matches_3D[i][k]] - 1;
      }
    }
    vtkNew<vtkExtractVOI> extractVOI;
    extractVOI->SetInputData(imageIn);
    extractVOI->SetVOI(VOI[0], VOI[1], VOI[2], VOI[3], VOI[4], VOI[5]);
    extractVOI->Update();
    vtkImageData *extracted = extractVOI->GetOutput();
    vtkSmartPointer<vtkCharArray> buffer = vtkSmartPointer<vtkCharArray>::New();
    if(vtkCommunicator::MarshalDataObject(extracted, buffer) == 0) {
      printErr("Marshalling failed!");
    };
    charArray3DBoundaries[matches_3D[i][0]].emplace_back(buffer);
    charArray3DBoundariesMetaData[matches_3D[i][0]].emplace_back(
      std::array<ttk::SimplexId, 3>{
        matches_3D[i][1], matches_3D[i][2], matches_3D[i][3]});
  }
  std::array<ttk::SimplexId, 3> sendMetadataArray3D;
  std::array<ttk::SimplexId, 3> recvMetaDataArray3D;
  std::vector<vtkSmartPointer<vtkCharArray>> charArray3DBoundariesReceived;
  std::vector<std::array<ttk::SimplexId, 3>>
    charArray3DBoundariesMetaDataReceived;
  for(int i = 0; i < ttk::MPIsize_; i++) {
    for(int j = 0; j < charArray3DBoundaries[i].size(); j++) {
      send_size = charArray3DBoundaries[i][j]->GetNumberOfTuples();
      MPI_Sendrecv(&send_size, 1, ttk::getMPIType(send_size), i, 0, &recv_size,
                   1, ttk::getMPIType(recv_size), i, 0, ttk::MPIcomm_,
                   MPI_STATUS_IGNORE);
      vtkSmartPointer<vtkCharArray> buffer
        = vtkSmartPointer<vtkCharArray>::New();
      buffer->SetNumberOfTuples(recv_size);
      buffer->SetNumberOfComponents(1);
      char *sendPointer
        = ttkUtils::GetPointer<char>(charArray3DBoundaries[i][j]);
      char *recvPointer = ttkUtils::GetPointer<char>(buffer);
      MPI_Sendrecv(sendPointer, send_size, MPI_CHAR, i, 0, recvPointer,
                   recv_size, MPI_CHAR, i, 0, ttk::MPIcomm_, MPI_STATUS_IGNORE);
      charArray3DBoundariesReceived.emplace_back(buffer);
      sendMetadataArray3D = charArray3DBoundariesMetaData[i][j];
      MPI_Sendrecv(&sendMetadataArray3D, 3, ttk::getMPIType(sendMetadata), i, 0,
                   &recvMetaDataArray3D, 3, ttk::getMPIType(recvMetaData), i, 0,
                   ttk::MPIcomm_, MPI_STATUS_IGNORE);
      charArray3DBoundariesMetaDataReceived.emplace_back(recvMetaDataArray3D);
    }
  }

  /*
  for (int i = 0; i < charArray3DBoundariesReceived.size(); i++){
    printMsg("Size of received charArray:
  "+std::to_string(charArray3DBoundariesReceived[i]->GetNumberOfTuples()));
  }*/
  imageOut->DeepCopy(imageIn);

  // Merge in the x direction (x_low and x_high)
  for(int dir = 0; dir < 2; dir++) {
    auto it = std::find(charArrayBoundariesMetaDataReceived.begin(),
                        charArrayBoundariesMetaDataReceived.end(), other(dir));
    if(it != charArrayBoundariesMetaDataReceived.end()) {
      vtkNew<vtkStructuredPoints> id;
      vtkNew<vtkImageData> mergedImage;
      if(vtkCommunicator::UnMarshalDataObject(
           charArrayBoundariesReceived[std::distance(
             charArrayBoundariesMetaDataReceived.begin(), it)],
           id)
         == 0) {
        printErr("UnMarshaling failed!");
      };
      this->MergeImageAppendAndSlice(imageIn, id, mergedImage, dir);
      imageOut->DeepCopy(mergedImage);
    }
  }
  // Merge in the y direction
  for(int dir = 2; dir < 4; dir++) {
    auto it = std::find(charArrayBoundariesMetaDataReceived.begin(),
                        charArrayBoundariesMetaDataReceived.end(), other(dir));
    vtkNew<vtkImageData> mergedImage;
    if(it != charArrayBoundariesMetaDataReceived.end()) {
      vtkNew<vtkStructuredPoints> id;
      if(vtkCommunicator::UnMarshalDataObject(
           charArrayBoundariesReceived[std::distance(
             charArrayBoundariesMetaDataReceived.begin(), it)],
           id)
         == 0) {
        printErr("UnMarshaling failed!");
      };
      mergedImage->DeepCopy(id);
    }
    for(int dir_2D = 0; dir_2D < 2; dir_2D++) {
      auto it_2D
        = std::find(charArray2DBoundariesMetaDataReceived.begin(),
                    charArray2DBoundariesMetaDataReceived.end(),
                    std::array<ttk::SimplexId, 2>{other(dir_2D), other(dir)});
      if(it_2D != charArray2DBoundariesMetaDataReceived.end()) {
        vtkNew<vtkStructuredPoints> id_2D;
        vtkNew<vtkImageData> aux;
        if(vtkCommunicator::UnMarshalDataObject(
             charArray2DBoundariesReceived[std::distance(
               charArray2DBoundariesMetaDataReceived.begin(), it_2D)],
             id_2D)
           == 0) {
          printErr("UnMarshaling failed!");
        };
        this->MergeImageAppendAndSlice(mergedImage, id_2D, aux, dir_2D);
        mergedImage->DeepCopy(aux);
      }
    }

    if(mergedImage->GetNumberOfPoints() > 0) {
      vtkNew<vtkImageData> aux;
      this->MergeImageAppendAndSlice(imageOut, mergedImage, aux, dir);
      imageOut->DeepCopy(aux);
    }
  }

  // Merge in the z direction

  printMsg("End of periodic ghost generation");
  return 1;
};

int ttkPeriodicGhostsGeneration::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  /*vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);
  vtkImageData* image = vtkImageData::SafeDownCast(input);*/
  this->MPIPeriodicGhostPipelinePreconditioning(inputVector, outputVector);
  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  /*ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;*/
  /*vtkImageData* image1 =
  vtkImageData::SafeDownCast(vtkDataSet::GetData(inputVector[0]));
  vtkSmartPointer<vtkImageAppend> imageAppend =
  vtkSmartPointer<vtkImageAppend>::New(); vtkNew<vtkImageData> input2; int*
  extent; double* origin; double origin2[3]; extent = image->GetExtent();
  origin = image->GetOrigin();*/
  // for (int i = 0; i < 2; i++){
  //  extent[i] += 5;
  //}
  /*origin2[0] = origin[0]+5;
  origin2[1] = origin[1];
  origin2[2] = origin[2];*/
  /*input2->DeepCopy(image);
  //input2->SetExtent(extent);
  //input2->SetOrigin(origin2);
  imageAppend->SetInputData(image);
  imageAppend->AddInputData(input2);
  imageAppend->SetAppendAxis(0);
  //imageAppend->PreserveExtentsOff();
  imageAppend->Update();*/

  /*vtkImageData* image = vtkImageData::GetData(inputVector[0]);
  vtkImageData *output = vtkImageData::GetData(outputVector);
  vtkNew<vtkImageAppend> imageAppend;
  imageAppend->SetInputData(image);
  imageAppend->AddInputData(image);
  imageAppend->SetAppendAxis(0);
  //mageAppend->PreserveExtentsOn();
  imageAppend->Update();
  //vtkInformation* outInfo = outputVector->GetInformationObject(0);
  //int* appendedExtent = imageAppend->GetOutput()->GetExtent();
  //outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
appendedExtent, 6); vtkNew<vtkImageData>outImage{}; int newExtent[6]; int*
extent = image->GetExtent(); newExtent[0] = extent[0]; newExtent[1] =
extent[1]+(extent[1]-extent[0]); newExtent[2] = extent[2]; newExtent[3] =
extent[3]; newExtent[4] = extent[4]; newExtent[5] = extent[5];
  outImage->SetExtent(newExtent);
  outImage->SetSpacing(image->GetSpacing());
  outImage->SetOrigin(image->GetOrigin());
  //outImage->SetNumber
  vtkNew<vtkCharArray> test;
  test->SetNumberOfTuples(outImage->GetNumberOfPoints());
  test->Fill(1);
  test->SetName("Test");
  //outInfo->Set(vtkStreamingDemandDrivenPipeline::UNRESTRICTED_UPDATE_EXTENT(),
1);
  //outInfo->Set(
  //  vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), newExtent, 6);
  outImage->GetPointData()->AddArray(test);
  output->ShallowCopy(imageAppend->GetOutput());
  outInfo->Set(vtkDataObject::SPACING(), outImage->GetSpacing(), 3);
  outInfo->Set(vtkDataObject::ORIGIN(), outImage->GetOrigin(), 3);
  //outInfo->Set(
  //  vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outImage->GetExtent(),
6); printMsg("inside :"+std::to_string(image->GetNumberOfPoints()));

  printMsg("outside after ghost :"+std::to_string(output->GetNumberOfPoints()));
/*  int* inWextent =
outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
  printMsg("inWextent: "+std::to_string(inWextent[0])+", "
  +std::to_string(inWextent[1])+", "
  +std::to_string(inWextent[2])+", "
  +std::to_string(inWextent[3])+", "
  +std::to_string(inWextent[4])+", "
  +std::to_string(inWextent[5]));
  printMsg("appendedExtent: "+std::to_string(appendedExtent[0])+", "
  +std::to_string(appendedExtent[1])+", "
  +std::to_string(appendedExtent[2])+", "
  +std::to_string(appendedExtent[3])+", "
  +std::to_string(appendedExtent[4])+", "
  +std::to_string(appendedExtent[5]));
  appendedExtent = image->GetExtent();
  printMsg("originExtent: "+std::to_string(appendedExtent[0])+", "
  +std::to_string(appendedExtent[1])+", "
  +std::to_string(appendedExtent[2])+", "
  +std::to_string(appendedExtent[3])+", "
  +std::to_string(appendedExtent[4])+", "
  +std::to_string(appendedExtent[5]));  */
  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  // vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);*/

  // make a SHALLOW copy of the input
  // outputDataSet->ShallowCopy(inputDataSet);

  // return success
  return 1;
};
