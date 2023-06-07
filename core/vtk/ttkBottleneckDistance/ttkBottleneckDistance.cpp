#include <ttkBottleneckDistance.h>
#include <ttkBottleneckDistanceUtils.h>
#include <ttkUtils.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkBottleneckDistance);

ttkBottleneckDistance::ttkBottleneckDistance() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

int ttkBottleneckDistance::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkBottleneckDistance::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

static int generateMatchings(vtkUnstructuredGrid *const outputCT3,
                             const ttk::DiagramType &diagram1,
                             const ttk::DiagramType &diagram2,
                             const std::vector<ttk::MatchingType> &matchings,
                             const std::array<double, 3> &distances,
                             const double globalDist,
                             const double spacing,
                             const bool isBottleneck,
                             const bool is2D0,
                             const bool is2D1) {

  vtkNew<vtkUnstructuredGrid> vtu{};

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(2 * matchings.size());
  vtu->SetPoints(points);

  vtkNew<vtkDoubleArray> costs{};
  costs->SetName("Cost");
  costs->SetNumberOfComponents(1);
  costs->SetNumberOfTuples(matchings.size());
  vtu->GetCellData()->AddArray(costs);

  vtkNew<vtkIntArray> matchingIds{};
  matchingIds->SetName("MatchingIdentifier");
  matchingIds->SetNumberOfComponents(1);
  matchingIds->SetNumberOfTuples(matchings.size());
  vtu->GetCellData()->AddArray(matchingIds);

  // Build matchings.
  for(size_t i = 0; i < matchings.size(); ++i) {
    const auto &t = matchings[i];
    const auto n1 = std::get<0>(t);
    const auto n2 = std::get<1>(t);

    const auto &pair0 = diagram1[n1];
    const auto &pair1 = diagram2[n2];

    const auto pairPoint = [](const ttk::PersistencePair &pair, const bool is2D,
                              const double zval) -> std::array<double, 3> {
      if(is2D) {
        return {pair.birth.sfValue, pair.death.sfValue, zval};
      } else {
        return {
          pair.death.coords[0], pair.death.coords[1], pair.death.coords[2]};
      }
    };

    const auto p0 = pairPoint(pair0, is2D0, -spacing / 2.0);
    points->SetPoint(2 * i + 0, p0.data());
    const auto p1 = pairPoint(pair1, is2D1, spacing / 2.0);
    points->SetPoint(2 * i + 1, p1.data());

    std::array<vtkIdType, 2> ids{
      2 * static_cast<vtkIdType>(i) + 0,
      2 * static_cast<vtkIdType>(i) + 1,
    };
    vtu->InsertNextCell(VTK_LINE, 2, ids.data());

    costs->SetTuple1(i, std::get<2>(t));
    matchingIds->SetTuple1(i, i);
  }

  // add distance results to output_matchings FieldData
  vtkNew<vtkDoubleArray> minSad{};
  minSad->SetName("MinSaddleCost");
  minSad->SetNumberOfTuples(1);
  minSad->SetTuple1(0, distances[0]);

  vtkNew<vtkDoubleArray> sadSad{};
  sadSad->SetName("SaddleSaddleCost");
  sadSad->SetNumberOfTuples(1);
  sadSad->SetTuple1(0, distances[1]);

  vtkNew<vtkDoubleArray> sadMax{};
  sadMax->SetName("SaddleMaxCost");
  sadMax->SetNumberOfTuples(1);
  sadMax->SetTuple1(0, distances[2]);

  vtkNew<vtkDoubleArray> wass{};
  wass->SetName(isBottleneck ? "BottleneckDistance" : "WassersteinDistance");
  wass->SetNumberOfTuples(1);
  wass->SetTuple1(0, globalDist);

  vtu->GetFieldData()->AddArray(minSad);
  vtu->GetFieldData()->AddArray(sadSad);
  vtu->GetFieldData()->AddArray(sadMax);
  vtu->GetFieldData()->AddArray(wass);

  outputCT3->ShallowCopy(vtu);

  return 1;
}

int ttkBottleneckDistance::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto outputDiagrams = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto outputMatchings = vtkUnstructuredGrid::GetData(outputVector, 1);

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0]);
  std::vector<vtkUnstructuredGrid *> inputDiags{};

  if(blocks == nullptr) {
    this->printErr("No input diagrams");
    return 0;
  }

  for(size_t i = 0; i < blocks->GetNumberOfBlocks(); ++i) {
    const auto diag = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
    if(diag != nullptr) {
      inputDiags.emplace_back(diag);
    }
  }

  if(inputDiags.size() < 2) {
    this->printErr("Less than two input diagrams");
    return 0;
  }
  if(inputDiags.size() > 2) {
    this->printWrn("More than two input diagrams: "
                   + std::to_string(inputDiags.size()));
  }

  const auto coords0 = vtkFloatArray::SafeDownCast(
    inputDiags[0]->GetPointData()->GetArray("Coordinates"));
  const auto coords1 = vtkFloatArray::SafeDownCast(
    inputDiags[1]->GetPointData()->GetArray("Coordinates"));

  const bool is2D0 = coords0 != nullptr;
  const bool is2D1 = coords1 != nullptr;

  // Call package
  int status = 0;

  ttk::DiagramType diagram0{}, diagram1{};
  std::vector<ttk::MatchingType> matchings{};

  status = VTUToDiagram(diagram0, inputDiags[0], *this);
  if(status < 0) {
    this->printErr("Could not extract diagram from first input data-set");
    return 0;
  }

  status = VTUToDiagram(diagram1, inputDiags[1], *this);
  if(status < 0) {
    this->printErr("Could not extract diagram from second input data-set");
    return 0;
  }

  // Exec.
  status = this->execute(diagram0, diagram1, matchings);
  if(status != 0) {
    this->printErr("Base layer failed with error status "
                   + std::to_string(status));
    return 0;
  }

  // Generate matchings
  if(this->UseOutputMatching) {
    status = generateMatchings(outputMatchings, diagram0, diagram1, matchings,
                               this->costs_, this->distance_, this->Spacing,
                               this->WassersteinMetric == "inf", is2D0, is2D1);

    if(status != 1) {
      this->printErr("Could not compute matchings");
      return 0;
    }
  }

  // Translate diagrams
  vtkNew<vtkUnstructuredGrid> vtu0{}, vtu1{};
  if(this->UseGeometricSpacing) {
    vtu0->ShallowCopy(inputDiags[0]);
    ResetDiagramPosition(vtu0, *this);
    TranslateDiagram(vtu0, {0, 0, -this->Spacing});
    vtu1->ShallowCopy(inputDiags[1]);
    ResetDiagramPosition(vtu1, *this);
    TranslateDiagram(vtu1, {0, 0, this->Spacing});
  } else {
    vtu0->ShallowCopy(inputDiags[0]);
    vtu1->ShallowCopy(inputDiags[1]);
  }

  // Add matchings infos on diagrams
  status = augmentDiagrams(matchings, vtu0, vtu1);
  if(status != 1) {
    this->printErr("Could not augment diagrams");
    return 0;
  }

  // Set output.
  outputDiagrams->SetNumberOfBlocks(2);
  outputDiagrams->SetBlock(0, vtu0);
  outputDiagrams->SetBlock(1, vtu1);

  return 1;
}
