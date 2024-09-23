#include <ttkBottleneckDistanceUtils.h>
#include <ttkMacros.h>
#include <ttkTrackingFromPersistenceDiagrams.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTrackingFromPersistenceDiagrams);

ttkTrackingFromPersistenceDiagrams::ttkTrackingFromPersistenceDiagrams() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkTrackingFromPersistenceDiagrams::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkTrackingFromPersistenceDiagrams::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkTrackingFromPersistenceDiagrams::buildMesh(
  const std::vector<ttk::trackingTuple> &trackings,
  const std::vector<std::vector<ttk::MatchingType>> &outputMatchings,
  const std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
  const bool useGeometricSpacing,
  const double spacing,
  const bool doPostProc,
  const std::vector<std::set<int>> &trackingTupleToMerged,
  vtkPoints *points,
  vtkUnstructuredGrid *persistenceDiagram,
  vtkDoubleArray *persistenceScalars,
  vtkDoubleArray *valueScalars,
  vtkIntArray *matchingIdScalars,
  vtkIntArray *lengthScalars,
  vtkIntArray *timeScalars,
  vtkIntArray *componentIds,
  vtkIntArray *pointTypeScalars,
  const ttk::Debug &dbg) {

  using ttk::CriticalType;
  int currentVertex = 0;
  for(size_t k = 0; k < trackings.size(); ++k) {
    const ttk::trackingTuple &tt = trackings[k];

    int const numStart = std::get<0>(tt);
    //     int numEnd = std::get<1>(tt);
    const std::vector<ttk::SimplexId> &chain = std::get<2>(tt);
    int const chainLength = chain.size();

    if(chainLength <= 1) {
      dbg.printErr("Got an unexpected 0-size chain.");
      return 0;
    }

    for(int c = 0; c < chainLength - 1; ++c) {
      const auto &matchings1 = outputMatchings[numStart + c];
      int const d1id = numStart + c;
      int const d2id = d1id + 1; // c % 2 == 0 ? d1id + 1 : d1id;
      const auto &diagram1 = inputPersistenceDiagrams[d1id];
      const auto &diagram2 = inputPersistenceDiagrams[d2id];

      // Insert segments
      vtkIdType ids[2];
      const auto n1 = chain[c];
      const auto n2 = chain[c + 1];

      // Search for right matching.
      double cost = 0.0;
      for(const auto &tup : matchings1) {
        auto d1id1 = (int)std::get<0>(tup);
        if(d1id1 == n1) {
          cost = std::get<2>(tup);
          break;
        }
      }

      const auto &pair0 = diagram1[n1];
      const auto &pair1 = diagram2[n2];

      double x1, y1, z1, x2, y2, z2;

      auto point1Type1 = pair0.birth.type;
      auto point1Type2 = pair0.death.type;
      auto point1Type = point1Type1 == CriticalType::Local_maximum
                            || point1Type2 == CriticalType::Local_maximum
                          ? CriticalType::Local_maximum
                        : point1Type1 == CriticalType::Local_minimum
                            || point1Type2 == CriticalType::Local_minimum
                          ? CriticalType::Local_minimum
                        : point1Type1 == CriticalType::Saddle2
                            || point1Type2 == CriticalType::Saddle2
                          ? CriticalType::Saddle2
                          : CriticalType::Saddle1;
      bool t11Min = point1Type1 == CriticalType::Local_minimum;
      bool t11Max = point1Type1 == CriticalType::Local_maximum;
      bool t12Min = point1Type2 == CriticalType::Local_minimum;
      bool t12Max = point1Type2 == CriticalType::Local_maximum;
      bool bothEx1 = (t11Min && t12Max) || (t11Max && t12Min);
      if(bothEx1) {
        x1 = t12Max ? pair0.death.coords[0] : pair0.birth.coords[0];
        y1 = t12Max ? pair0.death.coords[1] : pair0.birth.coords[1];
        z1 = t12Max ? pair0.death.coords[2] : pair0.birth.coords[2];
        if(useGeometricSpacing)
          z1 += spacing * (numStart + c);
      } else {
        x1 = t12Max   ? pair0.death.coords[0]
             : t11Min ? pair0.birth.coords[0]
                      : (pair0.death.coords[0] + pair0.birth.coords[0]) / 2.0;
        y1 = t12Max   ? pair0.death.coords[1]
             : t11Min ? pair0.birth.coords[1]
                      : (pair0.death.coords[1] + pair0.birth.coords[1]) / 2.0;
        z1 = t12Max   ? pair0.death.coords[2]
             : t11Min ? pair0.birth.coords[2]
                      : (pair0.death.coords[2] + pair0.birth.coords[2]) / 2.0;
        if(useGeometricSpacing)
          z1 += spacing * (numStart + c);
      }

      // Postproc component ids.
      int cid = k;
      bool hasMergedFirst = false;
      if(doPostProc) {
        const auto &connected = trackingTupleToMerged[k];
        if(!connected.empty()) {
          int const min = *(connected.begin());
          const ttk::trackingTuple &ttt = trackings[min];
          // int numStart2 = std::get<0>(ttt);
          int const numEnd2 = std::get<1>(ttt);
          if((numEnd2 > 0 && numStart + c > numEnd2 + 1) && min < (int)k) {
            dbg.printMsg("Switched " + std::to_string(k) + " for "
                         + std::to_string(min));
            cid = min;
            hasMergedFirst = numStart + c <= numEnd2 + 3;
          }

          if(hasMergedFirst) {
            // std::cout << "Has merged first " << std::endl;

            // Replace former first end of the segment with previous ending
            // segment.
            const std::vector<ttk::SimplexId> &chain3 = std::get<2>(ttt);
            const auto nn = chain3.back();
            const auto &diagramRematch = inputPersistenceDiagrams[numEnd2];
            const auto &pairN = diagramRematch[nn];

            point1Type1 = pairN.birth.type;
            point1Type2 = pairN.death.type;
            point1Type = point1Type1 == CriticalType::Local_maximum
                             || point1Type2 == CriticalType::Local_maximum
                           ? CriticalType::Local_maximum
                         : point1Type1 == CriticalType::Local_minimum
                             || point1Type2 == CriticalType::Local_minimum
                           ? CriticalType::Local_minimum
                         : point1Type1 == CriticalType::Saddle2
                             || point1Type2 == CriticalType::Saddle2
                           ? CriticalType::Saddle2
                           : CriticalType::Saddle1;
            t11Min = point1Type1 == CriticalType::Local_minimum;
            t11Max = point1Type1 == CriticalType::Local_maximum;
            t12Min = point1Type2 == CriticalType::Local_minimum;
            t12Max = point1Type2 == CriticalType::Local_maximum;
            bothEx1 = (t11Min && t12Max) || (t11Max && t12Min);
            // std::cout << "xyz " << x1 << ", " << y1 << ", " << z1 <<
            // std::endl;
            if(bothEx1) {
              x1 = t12Max ? pairN.death.coords[0] : pairN.birth.coords[0];
              y1 = t12Max ? pairN.death.coords[1] : pairN.birth.coords[1];
              z1 = t12Max ? pairN.death.coords[2] : pairN.birth.coords[2];
              if(useGeometricSpacing)
                z1 += spacing * (numStart + c);
            } else {
              x1 = t12Max ? pairN.death.coords[0]
                   : t11Min
                     ? pairN.birth.coords[0]
                     : (pairN.birth.coords[0] + pairN.death.coords[0]) / 2.0;
              y1 = t12Max ? pairN.death.coords[1]
                   : t11Min
                     ? pairN.birth.coords[1]
                     : (pairN.birth.coords[1] + pairN.death.coords[1]) / 2.0;
              z1 = t12Max ? pairN.death.coords[2]
                   : t11Min
                     ? pairN.birth.coords[2]
                     : (pairN.birth.coords[2] + pairN.death.coords[2]) / 2.0;
              if(useGeometricSpacing)
                z1 += spacing * (numStart + c);
            }
            // std::cout << "xyz " << x1 << ", " << y1 << ", " << z1 <<
            // std::endl;
          }
        }
      }

      points->InsertNextPoint(x1, y1, z1);
      ids[0] = 2 * currentVertex;
      pointTypeScalars->InsertTuple1(ids[0], (double)(int)point1Type);
      timeScalars->InsertTuple1(ids[0], (double)numStart + c);
      componentIds->InsertTuple1(ids[0], cid);

      const auto point2Type1 = pair1.birth.type;
      const auto point2Type2 = pair1.death.type;
      const auto point2Type = point2Type1 == CriticalType::Local_maximum
                                  || point2Type2 == CriticalType::Local_maximum
                                ? CriticalType::Local_maximum
                              : point2Type1 == CriticalType::Local_minimum
                                  || point2Type2 == CriticalType::Local_minimum
                                ? CriticalType::Local_minimum
                              : point2Type1 == CriticalType::Saddle2
                                  || point2Type2 == CriticalType::Saddle2
                                ? CriticalType::Saddle2
                                : CriticalType::Saddle1;
      bool const t21Ex = point2Type1 == CriticalType::Local_minimum
                         || point2Type1 == CriticalType::Local_maximum;
      bool const t22Ex = point2Type2 == CriticalType::Local_minimum
                         || point2Type2 == CriticalType::Local_maximum;
      bool const bothEx2 = t21Ex && t22Ex;
      if(bothEx2) {
        x2 = point2Type2 == CriticalType::Local_maximum ? pair1.death.coords[0]
                                                        : pair1.birth.coords[0];
        y2 = point2Type2 == CriticalType::Local_maximum ? pair1.death.coords[1]
                                                        : pair1.birth.coords[1];
        z2 = point2Type2 == CriticalType::Local_maximum ? pair1.death.coords[2]
                                                        : pair1.birth.coords[2];
        if(useGeometricSpacing)
          z2 += spacing * (numStart + c + 1);
      } else {
        x2 = t22Ex   ? pair1.death.coords[0]
             : t21Ex ? pair1.birth.coords[0]
                     : (pair1.birth.coords[0] + pair1.death.coords[0]) / 2;
        y2 = t22Ex   ? pair1.death.coords[1]
             : t21Ex ? pair1.birth.coords[1]
                     : (pair1.birth.coords[1] + pair1.death.coords[1]) / 2;
        z2 = t22Ex   ? pair1.death.coords[2]
             : t21Ex ? pair1.birth.coords[2]
                     : (pair1.birth.coords[2] + pair1.death.coords[2]) / 2;
        if(useGeometricSpacing)
          z2 += spacing * (numStart + c + 1);
      }
      points->InsertNextPoint(x2, y2, z2);
      ids[1] = 2 * currentVertex + 1;
      pointTypeScalars->InsertTuple1(ids[1], (double)(int)point2Type);
      timeScalars->InsertTuple1(ids[1], (double)numStart + c);

      // Postproc component ids.
      componentIds->InsertTuple1(ids[1], cid);

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);

      persistenceScalars->InsertTuple1(currentVertex, cost);
      valueScalars->InsertTuple1(
        currentVertex, (pair0.death.sfValue + pair1.death.sfValue) / 2);
      matchingIdScalars->InsertTuple1(currentVertex, currentVertex);
      lengthScalars->InsertTuple1(currentVertex, chainLength);

      currentVertex++;
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(valueScalars);
  persistenceDiagram->GetCellData()->AddArray(matchingIdScalars);
  persistenceDiagram->GetCellData()->AddArray(lengthScalars);
  persistenceDiagram->GetPointData()->AddArray(timeScalars);
  persistenceDiagram->GetPointData()->AddArray(componentIds);
  persistenceDiagram->GetPointData()->AddArray(pointTypeScalars);

  return 0;
}

int ttkTrackingFromPersistenceDiagrams::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // Unified bound fields
  vtkUnstructuredGrid *mesh = vtkUnstructuredGrid::GetData(outputVector, 0);

  // Number of input files
  int const numInputs = inputVector[0]->GetNumberOfInformationObjects();
  this->printMsg("Number of inputs: " + std::to_string(numInputs));

  // Get input data
  std::vector<vtkUnstructuredGrid *> inputVTUs(numInputs);
  for(int i = 0; i < numInputs; i++) {
    inputVTUs[i] = vtkUnstructuredGrid::GetData(inputVector[0], i);
  }

  std::vector<ttk::DiagramType> inputPersistenceDiagrams(numInputs);
  std::vector<vtkNew<vtkUnstructuredGrid>> outputDiags(2 * numInputs - 2);
  std::vector<std::vector<ttk::MatchingType>> outputMatchings(numInputs - 1);

  // Input parameters.
  double const spacing = Spacing;
  std::string const algorithm = DistanceAlgorithm;
  double const tolerance = Tolerance;
  std::string const wasserstein = WassersteinMetric;

  // Transform inputs into the right structure.
  for(int i = 0; i < numInputs; ++i) {
    VTUToDiagram(inputPersistenceDiagrams[i], inputVTUs[i], *this);
  }

  this->performMatchings(
    numInputs, inputPersistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, PX, PY, PZ, PS, PE // Coefficients
  );

  // Get back meshes.
  //  #pragma omp parallel for num_threads(ThreadNumber)
  for(int i = 0; i < numInputs - 1; ++i) {
    outputDiags[2 * i + 0]->ShallowCopy(inputVTUs[i]);
    outputDiags[2 * i + 1]->ShallowCopy(inputVTUs[i + 1]);

    int const status = augmentDiagrams(
      outputMatchings[i], outputDiags[2 * i + 0], outputDiags[2 * i + 1]);
    if(status < 0)
      return -2;
  }

  for(int i = 0; i < numInputs; ++i) {
    const auto &grid = outputDiags[i];
    if(!grid || !grid->GetCellData()
       || !grid->GetCellData()->GetArray("Persistence")) {
      this->printErr("Inputs are not persistence diagrams");
      return 0;
    }

    // Check if inputs have the same data type and the same number of points
    if(grid->GetCellData()->GetArray("Persistence")->GetDataType()
       != outputDiags[0]
            ->GetCellData()
            ->GetArray("Persistence")
            ->GetDataType()) {
      this->printErr("Inputs of different data types");
      return 0;
    }
  }

  for(int i = 0; i < numInputs - 1; ++i) {
    const auto &grid1 = outputDiags[i];
    const auto &grid2 = outputDiags[i + 1];
    if(i % 2 == 1 && i < numInputs - 1
       && grid1->GetCellData()->GetNumberOfTuples()
            != grid2->GetCellData()->GetNumberOfTuples()) {
      this->printErr("Inconsistent length or order of input diagrams.");
      return 0;
    }
  }

  vtkUnstructuredGrid *outputMesh = vtkUnstructuredGrid::SafeDownCast(mesh);

  vtkNew<vtkPoints> const points{};
  vtkNew<vtkUnstructuredGrid> const persistenceDiagram{};

  vtkNew<vtkDoubleArray> persistenceScalars{};
  vtkNew<vtkDoubleArray> valueScalars{};
  vtkNew<vtkIntArray> matchingIdScalars{};
  vtkNew<vtkIntArray> lengthScalars{};
  vtkNew<vtkIntArray> timeScalars{};
  vtkNew<vtkIntArray> componentIds{};
  vtkNew<vtkIntArray> pointTypeScalars{};
  persistenceScalars->SetName("Cost");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  // (+ vertex id)
  std::vector<ttk::trackingTuple>
    trackingsBase; // structure containing all trajectories
  this->performTracking(
    inputPersistenceDiagrams, outputMatchings, trackingsBase);

  std::vector<std::set<int>> trackingTupleToMerged(trackingsBase.size());
  if(DoPostProc)
    this->performPostProcess(inputPersistenceDiagrams, trackingsBase,
                             trackingTupleToMerged, PostProcThresh);

  // bool Is3D = true;
  bool const useGeometricSpacing = UseGeometricSpacing;
  // auto spacing = (float) Spacing;

  // std::vector<trackingTuple> trackings = *trackingsBase;
  // Row = iteration number
  // Col = (id in pd 2, arrival point id)

  // Build mesh.
  buildMesh(trackingsBase, outputMatchings, inputPersistenceDiagrams,
            useGeometricSpacing, spacing, DoPostProc, trackingTupleToMerged,
            points, persistenceDiagram, persistenceScalars, valueScalars,
            matchingIdScalars, lengthScalars, timeScalars, componentIds,
            pointTypeScalars, *this);

  outputMesh->ShallowCopy(persistenceDiagram);

  return 1;
}
