#include <ttkBottleneckDistance.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>

#include <random>

vtkStandardNewMacro(ttkBottleneckDistance);

ttkBottleneckDistance::ttkBottleneckDistance() {
  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(3);
}

int ttkBottleneckDistance::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(ttkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkBottleneckDistance::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int generatePersistenceDiagram(ttk::DiagramType &diagram, const int size) {

  int v0 = 1;
  int v1 = 2;

  std::default_random_engine generator1(1998985);
  std::default_random_engine generator2(8584584);
  std::uniform_real_distribution<> dis1(0.0, 1.0);
  std::uniform_real_distribution<> dis2(0.0, 1.0);

  for(int i = 1; i < size; ++i) {
    int r0 = 2; // (rand() % 3);
    float r1 = 0.001f + dis1(generator1);
    float r2 = 0.001f + dis2(generator2);

    // BLocalMin BSaddle1 BSaddle2 BLocalMax
    int pairType = r0; // (0/min, 1/saddle, 2/max)
    BNodeType nodeType1; //
    BNodeType nodeType2; //
    switch(pairType) {
      // case 0:
      //   nodeType1 = BLocalMin;
      //   nodeType2 = BSaddle1;
      //   break;
      // case 1:
      //   nodeType1 = BSaddle1;
      //   nodeType2 = BSaddle2;
      //   break;
      case 2:
      default:
        nodeType1 = BSaddle2;
        nodeType2 = BLocalMax;
        break;
        // default:
        //   nodeType1 = (BNodeType) -1;
        //   nodeType2 = (BNodeType) -1;
    }

    float x1 = 0.5f * r1;
    float y1 = x1;
    float z1 = 0.f;

    float x2 = x1;
    float y2 = x1 + 0.5f * r2; // x1 + rand(0.5)
    float z2 = 0.f; // 0

    const auto birth = x1;
    const auto death = y2;

    const auto pers = death - birth;

    diagram.push_back(std::make_tuple(v0, nodeType1, v1, nodeType2, pers,
                                      pairType, birth, x1, y1, z1, death, x2,
                                      y2, z2));

    v0++;
    v1++;
  }

  sort(diagram.begin(), diagram.end(),
       [](const ttk::PairTuple &a, const ttk::PairTuple &b) -> bool {
         return std::get<6>(a) < std::get<6>(b);
       });

  return 1;
}

// Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
template <typename dataType>
int ttkBottleneckDistance::augmentPersistenceDiagrams(
  const std::vector<diagramTuple> &diagram1,
  const std::vector<diagramTuple> &diagram2,
  const std::vector<matchingTuple> &matchings,
  vtkUnstructuredGrid *const CTPersistenceDiagram1,
  vtkUnstructuredGrid *const CTPersistenceDiagram2) {

  auto diagramSize1 = (BIdVertex)diagram1.size();
  auto diagramSize2 = (BIdVertex)diagram2.size();
  auto matchingsSize = (BIdVertex)matchings.size();

  vtkNew<vtkIntArray> matchingIdentifiers1{};
  matchingIdentifiers1->SetName("MatchingIdentifier");

  vtkNew<vtkIntArray> matchingIdentifiers2{};
  matchingIdentifiers2->SetName("MatchingIdentifier");

  if(matchingsSize > 0) {
    ttk::SimplexId ids[2];
    matchingIdentifiers1->SetNumberOfComponents(1);
    matchingIdentifiers2->SetNumberOfComponents(1);
    matchingIdentifiers1->SetNumberOfTuples(diagramSize1);
    matchingIdentifiers2->SetNumberOfTuples(diagramSize2);

    // Unaffected by default
    for(BIdVertex i = 0; i < diagramSize1; ++i)
      matchingIdentifiers1->InsertTuple1(i, -1);
    for(BIdVertex i = 0; i < diagramSize2; ++i)
      matchingIdentifiers2->InsertTuple1(i, -1);

    // Last cell = junction
    if(diagramSize1
       < CTPersistenceDiagram1->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers1->InsertTuple1(diagramSize1, -1);
      matchingIdentifiers1->InsertTuple1(diagramSize1 + 1, -1);
    }
    if(diagramSize2
       < CTPersistenceDiagram2->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers2->InsertTuple1(diagramSize2, -1);
      matchingIdentifiers2->InsertTuple1(diagramSize2 + 1, -1);
    }

    // Affect bottleneck matchings
    int pairingIndex = 0;
    for(BIdVertex i = 0; i < matchingsSize; ++i) {
      matchingTuple t = matchings.at((unsigned long)i);
      ids[0] = std::get<0>(t);
      ids[1] = std::get<1>(t);
      matchingIdentifiers1->InsertTuple1(ids[0], pairingIndex);
      matchingIdentifiers2->InsertTuple1(ids[1], pairingIndex);
      pairingIndex++;
    }

    CTPersistenceDiagram1->GetCellData()->AddArray(matchingIdentifiers1);
    CTPersistenceDiagram2->GetCellData()->AddArray(matchingIdentifiers2);
  }

  return 1;
}

int translateDiagram(vtkUnstructuredGrid *output,
                     vtkUnstructuredGrid *input,
                     const double spacing) {
  vtkNew<vtkTransform> tr{};
  tr->Translate(0, 0, spacing);
  vtkNew<vtkTransformFilter> trf{};
  trf->SetTransform(tr);
  trf->SetInputData(input);
  trf->Update();
  output->ShallowCopy(trf->GetOutputDataObject(0));

  return 1;
}

template <typename dataType>
int getMatchingMesh(vtkUnstructuredGrid *const outputCT3,
                    const ttk::DiagramType &diagram1,
                    const ttk::DiagramType &diagram2,
                    const std::vector<matchingTuple> &matchings,
                    const double spacing,
                    const bool is2D) {

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

    const auto pairPoint = [is2D](const ttk::PairTuple &pair,
                                  const double zval) -> std::array<double, 3> {
      if(is2D) {
        return {std::get<6>(pair), std::get<10>(pair), zval};
      } else {
        return {std::get<11>(pair), std::get<12>(pair), std::get<13>(pair)};
      }
    };

    const auto p0 = pairPoint(pair0, spacing / 2.0);
    points->SetPoint(2 * i + 0, p0.data());
    const auto p1 = pairPoint(pair1, -spacing / 2.0);
    points->SetPoint(2 * i + 1, p1.data());

    std::array<vtkIdType, 2> ids{
      2 * static_cast<vtkIdType>(i) + 0,
      2 * static_cast<vtkIdType>(i) + 1,
    };
    vtu->InsertNextCell(VTK_LINE, 2, ids.data());

    costs->SetTuple1(i, std::get<2>(t));
    matchingIds->SetTuple1(i, i);
  }

  outputCT3->ShallowCopy(vtu);

  return 1;
}

int ttkBottleneckDistance::doBenchmark() {
  ttk::DiagramType CTDiagram1{}, CTDiagram2{};

  int status = 0;
  status = generatePersistenceDiagram(CTDiagram1, this->BenchmarkSize);
  if(status < 0)
    return status;
  status = generatePersistenceDiagram(CTDiagram2, 4 * this->BenchmarkSize);
  if(status < 0)
    return status;

  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);
  this->setCTDiagram1(&CTDiagram1);
  this->setCTDiagram2(&CTDiagram2);

  this->setWasserstein(this->WassersteinMetric);
  this->setAlgorithm(this->DistanceAlgorithm);
  this->setPVAlgorithm(this->PVAlgorithm);

  // Empty matchings.
  ttk::DiagramType matchings{};
  this->setOutputMatchings(&matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  status = this->execute<double>(usePersistenceMetric);

  return status;
}

int ttkBottleneckDistance::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  if(this->BenchmarkSize > 0) {
    return doBenchmark();
  }

  auto outputCT1 = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputCT2 = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputCT3 = vtkUnstructuredGrid::GetData(outputVector, 2);

  // Wrap
  // this->setWrapper(this);
  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);

  auto inputDiag0 = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto inputDiag1 = vtkUnstructuredGrid::GetData(inputVector[1]);

  if(!inputDiag0 || !inputDiag1 || !outputCT3) {
    this->printErr("Input grids should be non-NULL");
    return 0;
  }

  auto birthScalars1 = vtkDoubleArray::SafeDownCast(
    inputDiag0->GetPointData()->GetArray("Birth"));
  auto deathScalars1 = vtkDoubleArray::SafeDownCast(
    inputDiag0->GetPointData()->GetArray("Death"));
  auto birthScalars2 = vtkDoubleArray::SafeDownCast(
    inputDiag1->GetPointData()->GetArray("Birth"));
  auto deathScalars2 = vtkDoubleArray::SafeDownCast(
    inputDiag1->GetPointData()->GetArray("Death"));
  bool is2D1 = !deathScalars1 && !birthScalars1;
  bool is2D2 = !deathScalars2 && !birthScalars2;
  if(is2D1 != is2D2) {
    this->printErr("Diagrams should not be embedded");
    return 0;
  }
  bool is2D = is2D1;

  // Call package
  int status = 0;

  ttk::DiagramType diagram0{}, diagram1{};

  status = VTUToDiagram(diagram0, inputDiag0, *this);
  if(status < 0) {
    this->printErr("Could not extract diagram from first input data-set");
    return 0;
  }

  status = VTUToDiagram(diagram1, inputDiag1, *this);
  if(status < 0) {
    this->printErr("Could not extract diagram from second input data-set");
    return 0;
  }

  this->setCTDiagram1(&diagram0);
  this->setCTDiagram2(&diagram1);

  this->setWasserstein(this->WassersteinMetric);
  this->setAlgorithm(this->DistanceAlgorithm);
  this->setPVAlgorithm(this->PVAlgorithm);

  // Empty matchings.
  std::vector<matchingTuple> matchings;
  this->setOutputMatchings(&matchings);

  // Exec.
  status = this->execute<double>(this->UsePersistenceMetric);
  if(status != 0) {
    this->printErr("Base layer failed with error status "
                   + std::to_string(status));
    return 0;
  }

  // Apply results to outputs 0 and 1.
  status = augmentPersistenceDiagrams<double>(
    diagram0, diagram1, matchings, inputDiag0, inputDiag1);
  if(status != 1) {
    this->printErr("Could not augment diagrams");
    return 0;
  }

  // Apply results to output 2.
  if(this->UseOutputMatching) {
    status = getMatchingMesh<double>(
      outputCT3, diagram0, diagram1, matchings, this->Spacing, is2D);

    if(status != 1) {
      this->printErr("Could not compute matchings");
      return 0;
    }
  }

  // Set output.
  if(this->UseGeometricSpacing) {
    translateDiagram(outputCT1, inputDiag0, this->Spacing / 2.0);
    translateDiagram(outputCT2, inputDiag1, -this->Spacing / 2.0);
  } else {
    outputCT1->ShallowCopy(inputDiag0);
    outputCT2->ShallowCopy(inputDiag1);
  }

  return 1;
}
