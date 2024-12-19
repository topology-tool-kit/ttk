#include <ttkContourTree.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

// VTK includes
#include <vtkConnectivityFilter.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyData.h>
#include <vtkThreshold.h>
#include <vtkVersionMacros.h>

vtkStandardNewMacro(ttkContourTree);

ttkContourTree::ttkContourTree() {
  this->setDebugMsgPrefix("ContourTree");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
  this->params_.treeType = ttk::ftm::TreeType::Contour;
}

int ttkContourTree::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkContourTree::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkContourTree::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputSkeletonNodes = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputSkeletonArcs = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 2);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Error: input pointer is NULL.");
    return 0;
  }

  if(!input->GetNumberOfPoints()) {
    this->printErr("Error: input has no point.");
    return 0;
  }

  if(!outputSkeletonNodes || !outputSkeletonArcs || !outputSegmentation) {
    this->printErr("Error: output pointer is NULL.");
    return 0;
  }
#endif

  // Arrays

  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray)
    return 0;

  // Connected components
  if(input->IsA("vtkUnstructuredGrid")) {
    // This data set may have several connected components,
    // we need to apply the FTM Tree for each one of these components
    // We then reconstruct the global tree using an offset mechanism
    auto inputWithId = vtkSmartPointer<vtkUnstructuredGrid>::New();
    inputWithId->ShallowCopy(input);
    identify(inputWithId);

    vtkNew<vtkConnectivityFilter> connectivity{};
    connectivity->SetInputData(inputWithId);
    connectivity->SetExtractionModeToAllRegions();
    connectivity->ColorRegionsOn();
    connectivity->Update();

    nbCC_ = connectivity->GetOutput()
              ->GetCellData()
              ->GetArray("RegionId")
              ->GetRange()[1]
            + 1;
    connected_components_.resize(nbCC_);

    if(nbCC_ > 1) {
      // Warning, in case of several connected components, the ids seen by
      // the base code will not be consistent with those of the original
      // mesh
      for(int cc = 0; cc < nbCC_; cc++) {
        vtkNew<vtkThreshold> threshold{};
        threshold->SetInputConnection(connectivity->GetOutputPort());
        threshold->SetInputArrayToProcess(
          0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
#if VTK_VERSION_NUMBER < VTK_VERSION_CHECK(9, 2, 0)
        threshold->ThresholdBetween(cc, cc);
#else
        threshold->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        threshold->SetLowerThreshold(cc);
        threshold->SetUpperThreshold(cc);
#endif
        threshold->Update();
        connected_components_[cc] = vtkSmartPointer<vtkUnstructuredGrid>::New();
        connected_components_[cc]->ShallowCopy(threshold->GetOutput());
      }
    } else {
      connected_components_[0] = inputWithId;
    }
  } else if(input->IsA("vtkPolyData")) {
    // NOTE: CC check should not be implemented on a per vtk module layer.
    nbCC_ = 1;
    connected_components_.resize(nbCC_);
    connected_components_[0] = vtkSmartPointer<vtkPolyData>::New();
    connected_components_[0]->ShallowCopy(input);
    identify(connected_components_[0]);
  } else {
    nbCC_ = 1;
    connected_components_.resize(nbCC_);
    connected_components_[0] = vtkSmartPointer<vtkImageData>::New();
    connected_components_[0]->ShallowCopy(input);
    identify(connected_components_[0]);
  }

  // now proceed for each triangulation obtained.

  if(preconditionTriangulation() == 0) {
    this->printErr("Error : wrong triangulation.");
    return 0;
  }

  // Fill the vector of scalar/offset, cut the array in pieces if needed
  if(getScalars() == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("Error : wrong input scalars.");
    return 0;
#endif
  }
  getOffsets();

  this->printMsg("Launching on field "
                 + std::string{inputScalars_[0]->GetName()});

  ttk::ftm::idNode acc_nbNodes = 0;

  // Build tree
  for(int cc = 0; cc < nbCC_; cc++) {
    ftmTree_[cc].tree.setVertexScalars(
      ttkUtils::GetVoidPointer(inputScalars_[cc]));
    ftmTree_[cc].tree.setVertexSoSoffsets(offsets_[cc].data());
    ftmTree_[cc].tree.setTreeType(ttk::ftm::TreeType::Contour);
    ftmTree_[cc].tree.setSegmentation(GetWithSegmentation());
    ftmTree_[cc].tree.setNormalizeIds(GetWithNormalize());

    ttkVtkTemplateMacro(inputArray->GetDataType(),
                        triangulation_[cc]->getType(),
                        (ftmTree_[cc].tree.build<VTK_TT, TTK_TT>(
                          (TTK_TT *)triangulation_[cc]->getData())));

    ftmTree_[cc].offset = acc_nbNodes;
    acc_nbNodes += ftmTree_[cc]
                     .tree.getTree(ttk::ftm::TreeType::Contour)
                     ->getNumberOfNodes();
  }

  UpdateProgress(0.50);

  // Construct output
  if(getSkeletonNodes(outputSkeletonNodes) == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("Error : wrong properties on skeleton nodes.");
    return 0;
#endif
  }

  if(getSkeletonArcs(outputSkeletonArcs) == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("Error : wrong properties on skeleton arcs.");
    return 0;
#endif
  }

  if(GetWithSegmentation()) {
    outputSegmentation->ShallowCopy(input);
    if(getSegmentation(outputSegmentation) == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
      this->printErr("Error : wrong properties on segmentation.");
      return 0;
#endif
    }
  }

  UpdateProgress(1);

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
  printCSVStats();
#endif

  return 1;
}

int ttkContourTree::getOffsets() {
  // should be called after getScalars for inputScalars_ needs to be filled

  offsets_.resize(nbCC_);
  for(int cc = 0; cc < nbCC_; cc++) {
    const auto offsets
      = this->GetOrderArray(connected_components_[cc], 0, triangulation_[cc],
                            false, 1, ForceInputOffsetScalarField);

    offsets_[cc].resize(connected_components_[cc]->GetNumberOfPoints());

    for(size_t i = 0; i < offsets_[cc].size(); i++) {
      offsets_[cc][i] = offsets->GetTuple1(i);
    }

#ifndef TTK_ENABLE_KAMIKAZE
    if(offsets_[cc].empty()) {
      this->printMsg(
        {"Error : wrong input offset scalar field for ", std::to_string(cc)},
        ttk::debug::Priority::ERROR);
      return 0;
    }
#endif
  }

  return 1;
}

int ttkContourTree::getScalars() {
  inputScalars_.resize(nbCC_);
  for(int cc = 0; cc < nbCC_; cc++) {
    inputScalars_[cc]
      = this->GetInputArrayToProcess(0, connected_components_[cc]);

#ifndef TTK_ENABLE_KAMIKAZE
    if(!inputScalars_[cc]) {
      this->printMsg({"Error : input scalar ", std::to_string(cc),
                      " field pointer is null."},
                     ttk::debug::Priority::ERROR);
      return 0;
    }
#endif
  }

  return 1;
}

int ttkContourTree::preconditionTriangulation() {
  triangulation_.resize(nbCC_);
  ftmTree_ = std::vector<ttk::ftm::LocalFTM>(nbCC_);

  for(int cc = 0; cc < nbCC_; cc++) {
    triangulation_[cc]
      = ttkAlgorithm::GetTriangulation(connected_components_[cc]);
    if(!triangulation_[cc]) {
      this->printErr("Error : ttkTriangulation::getTriangulation() is null.");
      return 0;
    }

    ftmTree_[cc].tree.setDebugLevel(debugLevel_);
    ftmTree_[cc].tree.setThreadNumber(threadNumber_);
    ftmTree_[cc].tree.preconditionTriangulation(triangulation_[cc]);

#ifndef TTK_ENABLE_KAMIKAZE
    if(triangulation_[cc]->isEmpty()) {
      this->printMsg({"Error : ttkTriangulation on connected component",
                      std::to_string(cc), " allocation problem."},
                     ttk::debug::Priority::ERROR);
      return 0;
    }
#endif
  }
  return 1;
}

// protected
