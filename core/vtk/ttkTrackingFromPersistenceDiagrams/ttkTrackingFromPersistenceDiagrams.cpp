#include                  <ttkTrackingFromPersistenceDiagrams.h>

vtkStandardNewMacro(ttkTrackingFromPersistenceDiagrams)

ttkTrackingFromPersistenceDiagrams::ttkTrackingFromPersistenceDiagrams()
{
  outputMesh_ = nullptr;
  UseAllCores = false;

  DistanceAlgorithm = "ttk";
  PVAlgorithm = -1;
  Alpha = 1.0;
  Tolerance = 1.0;
  DoPostProc = false;
  PostProcThresh = 0.0;
  PX = 1;
  PY = 1;
  PZ = 1;
  PE = 1;
  PS = 1;

  WassersteinMetric = "1";
  UseGeometricSpacing = false;
  Is3D = false;
  Spacing = 1.0;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

ttkTrackingFromPersistenceDiagrams::~ttkTrackingFromPersistenceDiagrams()
{
  if (outputMesh_)
    outputMesh_->Delete();
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkTrackingFromPersistenceDiagrams::needsToAbort()
{
  return GetAbortExecute() > 0;
}

// transmit progress status -- to copy paste in other wrappers
int ttkTrackingFromPersistenceDiagrams::updateProgress(const float &progress)
{
  {
    std::stringstream msg;
    msg << "[ttkTrackingFromPersistenceDiagrams] " << progress*100
      << "% processed...." << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkTrackingFromPersistenceDiagrams::doIt(
  vtkDataSet **input,
  vtkUnstructuredGrid *mesh,
  int numInputs)
{
  using dataType = double;

  auto inputPersistenceDiagrams =
    new std::vector<std::vector<diagramTuple>*>((unsigned long) numInputs);
  auto outputMatchings =
    new std::vector<std::vector<matchingTuple>*>((unsigned long) numInputs - 1);
  auto outputPersistenceDiagrams =
    new std::vector<vtkUnstructuredGrid*>((unsigned long) 2 * numInputs - 2);

  // std::vector<vtkUnstructuredGrid*> input;
  // for (int i = 0; i < numInputs; ++i) {
  //   vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(inputs[i]);
  //   input.push_back(grid);
  // }
  // auto cmp = [](vtkUnstructuredGrid* a,
  //               vtkUnstructuredGrid* b)
  // { return a->GetMTime() > b->GetMTime(); };
  // std::sort(input.begin(), input.end(), cmp);

  // Input parameters.
  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
//  int threadNumber = ThreadNumber;
//  int benchmark = 0;
  double alpha = Alpha;
  double tolerance = Tolerance;
//  bool useSpacing = UseGeometricSpacing;
  bool is3D = Is3D;
  std::string wasserstein = WassersteinMetric;
//  bool useMatchingMesh = false;
//  bool persistence = false;

  // Transform inputs into the right structure.
  for (int i = 0; i < numInputs; ++i) {
    // TODO do it in another, normalized loop.
    vtkSmartPointer<ttkBottleneckDistance> bottleneckModule = vtkSmartPointer<ttkBottleneckDistance>::New();
    vtkUnstructuredGrid* grid1 = vtkUnstructuredGrid::New();
    grid1->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i]));

    std::vector<diagramTuple>* CTDiagram1 = new std::vector<diagramTuple>();
    bottleneckModule->getPersistenceDiagram(CTDiagram1, grid1, spacing, 0);
    inputPersistenceDiagrams->at((unsigned long) i) = CTDiagram1;
  }

  tracking_.performMatchings(
    numInputs,
    inputPersistenceDiagrams,
    outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein,
    tolerance,
    is3D,
    alpha, // Blending
    PX, PY, PZ, PS, PE, // Coefficients
    this // Wrapper for accessing threadNumber
  );

  // Get back meshes.
//  #pragma omp parallel for num_threads(ThreadNumber)
  for (int i = 0; i < numInputs - 1; ++i)
  {
    vtkUnstructuredGrid* grid1 = vtkUnstructuredGrid::New();
    grid1->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i]));
    vtkUnstructuredGrid* grid2 = vtkUnstructuredGrid::New();
    grid2->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i + 1]));
    // std::cout << "i - " << i << std::endl;
    // std::cout << grid1->GetFieldData()->GetArray(0)->GetTuple1(0) << std::endl;
    // std::cout << grid2->GetFieldData()->GetArray(0)->GetTuple1(0) << std::endl;

    vtkSmartPointer<ttkBottleneckDistance> bottleneckModule = vtkSmartPointer<ttkBottleneckDistance>::New();
//    bottleneckModule->AddInputData(0, grid1);
//    bottleneckModule->AddInputData(1, grid2);
//    bottleneckModule->SetDistanceAlgorithm(algorithm);
//    bottleneckModule->SetPVAlgorithm(PVAlgorithm);
//    bottleneckModule->SetWassersteinMetric(wasserstein);
//    bottleneckModule->SetAlpha(alpha);
//    bottleneckModule->SetTolerance(tolerance);
//    bottleneckModule->SetPX(PX);
//    bottleneckModule->SetPY(PY);
//    bottleneckModule->SetPZ(PZ);
//    bottleneckModule->SetPE(PE);
//    bottleneckModule->SetPS(PS);

    // Unused.
//    bottleneckModule->SetUsePersistenceMetric(persistence);
//    bottleneckModule->SetUseOutputMatching(useMatchingMesh);
//    bottleneckModule->SetBenchmarkSize(benchmark);

//    bottleneckModule->SetUseGeometricSpacing(useSpacing);
//    bottleneckModule->SetIs3D(is3D);
//    bottleneckModule->SetSpacing(spacing);

//    bottleneckModule->SetThreadNumber(threadNumber);
//    bottleneckModule->Update();

    std::vector<diagramTuple>* CTDiagram1 = inputPersistenceDiagrams->at((unsigned long) i);
    std::vector<diagramTuple>* CTDiagram2 = inputPersistenceDiagrams->at((unsigned long) i + 1);
    std::vector<matchingTuple>* Matchings12 = outputMatchings->at((unsigned long) i);

    vtkUnstructuredGrid *CTPersistenceDiagram1_ = vtkUnstructuredGrid::SafeDownCast(grid1);
    vtkUnstructuredGrid *CTPersistenceDiagram2_ = vtkUnstructuredGrid::SafeDownCast(grid2);

    int status = bottleneckModule->augmentPersistenceDiagrams<dataType>(
      CTDiagram1, CTDiagram2, Matchings12, CTPersistenceDiagram1_, CTPersistenceDiagram2_);
    if (status < 0) return -2;

//    vtkUnstructuredGrid* output1 = vtkUnstructuredGrid::SafeDownCast(bottleneckModule->GetOutput(0));
//    vtkUnstructuredGrid* output2 = vtkUnstructuredGrid::SafeDownCast(bottleneckModule->GetOutput(1));

    // TODO [] copy input persistence diags into var & convert formats
    // TODO augmentPersistenceDiagrams
    // TODO getMatchingMesh
    // TODO create outputs

    // Question: tracking from fields no basecode? modular?

    outputPersistenceDiagrams->at((unsigned long) 2 * i) = vtkUnstructuredGrid::New();
    outputPersistenceDiagrams->at((unsigned long) 2 * i + 1) = vtkUnstructuredGrid::New();
    outputPersistenceDiagrams->at((unsigned long) 2 * i)->ShallowCopy(CTPersistenceDiagram1_);
    outputPersistenceDiagrams->at((unsigned long) 2 * i + 1)->ShallowCopy(CTPersistenceDiagram2_);
//    outputPersistenceDiagrams->at((unsigned long) 2 * i)->ShallowCopy(output1);
//    outputPersistenceDiagrams->at((unsigned long) 2 * i + 1)->ShallowCopy(output2);
  }

  int numPersistenceDiagramsInput = (int) outputPersistenceDiagrams->size(); // numInputs;

  for (int i = 0; i < numPersistenceDiagramsInput; ++i)
  {
    vtkUnstructuredGrid* grid = outputPersistenceDiagrams->at((unsigned long) i);
    if (!grid || !grid->GetCellData() ||
        !grid->GetCellData()->GetArray("Persistence"))
    {
      std::stringstream msg;
      msg << "[ttkTrackingFromPersistenceDiagrams] Inputs are not persistence diagrams." << std::endl;
      dMsg(cerr, msg.str(), fatalMsg);
      return -1;
    }

    // Check if inputs have the same data type and the same number of points
    if (grid->GetCellData()->GetArray("Persistence")->GetDataType()
        != outputPersistenceDiagrams->at(0)->GetCellData()->GetArray("Persistence")->GetDataType())
    {
      std::stringstream msg;
      msg << "[ttkTrackingFromPersistenceDiagrams] Inputs of different data types." << std::endl;
      dMsg(cerr, msg.str(), fatalMsg);
      return -3;
    }
  }

  for (int i = 0; i < numPersistenceDiagramsInput - 1; ++i)
  {
    vtkUnstructuredGrid* grid1 = outputPersistenceDiagrams->at((unsigned long) i);
    vtkUnstructuredGrid* grid2 = outputPersistenceDiagrams->at((unsigned long) i + 1);
    if (i % 2 == 1 && i < numPersistenceDiagramsInput - 1 &&
        grid1->GetCellData()->GetNumberOfTuples() !=
        grid2->GetCellData()->GetNumberOfTuples())
    {
      std::stringstream msg;
      msg << "[ttkTrackingFromPersistenceDiagrams] Inconsistent length or order of input diagrams." << std::endl;
      dMsg(cerr, msg.str(), fatalMsg);
      return -2;
    }
  }

  outputMesh_ = vtkUnstructuredGrid::New();
  vtkUnstructuredGrid *outputMesh = vtkUnstructuredGrid::SafeDownCast(mesh);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> persistenceScalars = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> valueScalars = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> matchingIdScalars = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> timeScalars = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> componentIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> pointTypeScalars = vtkSmartPointer<vtkIntArray>::New();
  persistenceScalars->SetName("Cost");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("NodeType");

  // (+ vertex id)

  // TODO set prototype.
  auto trackingsBase = new std::vector<trackingTuple>(); // trackings;
  tracking_.performTracking(
    inputPersistenceDiagrams, outputMatchings,
    trackingsBase);

  auto trackingTupleToMerged = new std::vector<std::set<int>>(trackingsBase->size(), std::set<int>());
  if (DoPostProc)
    tracking_.performPostProcess(
      inputPersistenceDiagrams,
      trackingsBase,
      trackingTupleToMerged,
      PostProcThresh);

  // bool Is3D = true;
  bool useGeometricSpacing = UseGeometricSpacing;
//  auto spacing = (float) Spacing;

//  std::vector<std::vector<matchingTuple>> allMatchings;
//  auto *allDiagrams = new std::vector<std::vector<diagramTuple>>();
  std::vector<trackingTuple> trackings = *trackingsBase;
  // Row = iteration number
  // Col = (id in pd 2, arrival point id)

//  allDiagrams->resize((unsigned long) numPersistenceDiagramsInput);
//
//  for (int in = 0; in < numPersistenceDiagramsInput - 1; ++in)
//  {
//    auto diag1 = outputPersistenceDiagrams->at((unsigned long) in);
//    auto diag2 = outputPersistenceDiagrams->at((unsigned long) in + 1);
//    // auto diag3 = outputPersistenceDiagrams->at(in + 2);
//
//    auto *diagram1 = new std::vector<diagramTuple>();
//    auto *diagram2 = new std::vector<diagramTuple>();
//    // auto *diagram3 = new std::vector<diagramTuple>();
//
//
//    // Do get diagrams.
//    for (int id = 0; id < 2; ++id)
//    {
//      std::vector<diagramTuple> *currentDiag = id == 0 ? diagram1 : diagram2;
//      vtkUnstructuredGrid *inputDiag = outputPersistenceDiagrams->at((unsigned long) in + id);
//      vtkIntArray* vertexIdentifierScalars = vtkIntArray::SafeDownCast(inputDiag-> GetPointData()->GetArray("VertexIdentifier"));
//      vtkIntArray* nodeTypeScalars = vtkIntArray::SafeDownCast(inputDiag->GetPointData()->GetArray("NodeType"));
//      vtkIntArray* pairIdentifierScalars = vtkIntArray::SafeDownCast(inputDiag->GetCellData()->GetArray("PairIdentifier"));
//      vtkIntArray* extremumIndexScalars = vtkIntArray::SafeDownCast(inputDiag->GetCellData()->GetArray("PairType"));
//      vtkDoubleArray* currentPersistenceScalars = vtkDoubleArray::SafeDownCast(inputDiag->GetCellData()->GetArray("Persistence"));
//      vtkDoubleArray* birthScalars = vtkDoubleArray::SafeDownCast(inputDiag->GetPointData()->GetArray("Birth"));
//      vtkDoubleArray* deathScalars = vtkDoubleArray::SafeDownCast(inputDiag->GetPointData()->GetArray("Death"));
//      vtkIntArray* matchingIds = vtkIntArray::SafeDownCast(inputDiag->GetCellData()->GetArray("MatchingIdentifier"));
//
//      vtkPoints* currentPoints = inputDiag->GetPoints();
//      auto pairingsSize = (int) pairIdentifierScalars->GetNumberOfTuples();
//      // auto s = (float) 0.0;
//
//      if (!deathScalars != !birthScalars) {
//        std::stringstream msg;
//        msg << "[ttkTrackingFromPersistenceDiagrams] Death and birth criterion failing." << std::endl;
//        dMsg(std::cout, msg.str(), fatalMsg);
//        return -2;
//      }
//      // Is3D = !(!deathScalars && !birthScalars);
//      // if (!Is3D) s = spacing;
//
//      if (pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars || !nodeTypeScalars
//          || !currentPersistenceScalars || !extremumIndexScalars || !currentPoints || !matchingIds)
//      {
//        std::stringstream msg;
//        msg << "[ttkTrackingFromPersistenceDiagrams] Some arrays are absent." << std::endl;
//        dMsg(std::cout, msg.str(), fatalMsg);
//        return -2;
//      }
//
//      currentDiag->resize((unsigned long) pairingsSize);
//
//      int nbNonCompact = 0;
//
//      int realIndex = 0;
//      for (int i = 0; i < pairingsSize; ++i) {
//        int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
//        int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
//        int nodeType1 = nodeTypeScalars->GetValue(2 * i);
//        int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);
//
//        int pairIdentifier = pairIdentifierScalars->GetValue(i);
//        int pairType = extremumIndexScalars->GetValue(i);
//        double currentPersistence = currentPersistenceScalars->GetValue(i);
//
//        vtkCell *c = inputDiag->GetCell(i);
//        vtkIdList *l = c->GetPointIds();
//        int nbPInCell = (int) l->GetNumberOfIds();
//        if (nbPInCell != 2) std::cout << "Bad." << std::endl;
//
//        int index1 = 2 * i;
//        int index2 = index1 + 1;
//
//        index1 = (int) l->GetId(0);
//        index2 = (int) l->GetId(1);
//
//        double* coords1 = currentPoints->GetPoint(index1);
//        auto x1 = (float) coords1[0];
//        auto y1 = (float) coords1[1];
//        auto z1 = (float) coords1[2];
//
//        double* coords2 = currentPoints->GetPoint(index2);
//        auto x2 = (float) coords2[0];
//        auto y2 = (float) coords2[1];
//        auto z2 = (float) coords2[2];
//
//        double value1 = (!birthScalars) ? (double) x1 : birthScalars->GetValue(2 * i);
//        double value2 = (!deathScalars) ? (double) y2 : deathScalars->GetValue(2 * i + 1);
//
//        // if (pairIdentifier != -1 && pairIdentifier < pairingsSize)
//        currentDiag->at((unsigned long) realIndex) = std::make_tuple(
//          vertexId1, (BNodeType) nodeType1,
//          vertexId2, (BNodeType) nodeType2,
//          currentPersistence,
//          pairType,
//          value1, x1, y1, z1,
//          value2, x2, y2, z2);
//
//        realIndex++;
//
//        if (pairIdentifier >= pairingsSize) {
//          nbNonCompact++;
//          if (nbNonCompact == 0) {
//            std::stringstream msg;
//            msg << "[TTKBottleneckDistance] Diagram pair identifiers "
//                << "must be compact (not exceed the diagram size). "
//                << std::endl;
//            dMsg(std::cout, msg.str(), timeMsg);
//          }
//        }
//      }
//
//      if (nbNonCompact > 0) {
//        std::stringstream msg;
//        msg << "[TTKBottleneckDistance] Missed "
//          << nbNonCompact << " pairs due to non-compactness."
//          << std::endl;
//        dMsg(std::cout, msg.str(), timeMsg);
//      }
//    }
//
//    allDiagrams->at((unsigned long) in) = *diagram1;
//    allDiagrams->at((unsigned long) in + 1) = *diagram2;
//    // allDiagrams->at((unsigned long) in + 2) = *diagram3;
//
//    // Rebuild Matchings.
//    std::vector<matchingTuple> matchings1;
//    std::vector<matchingTuple> matchings2;
//    vtkIntArray* matchingIds1 = vtkIntArray::SafeDownCast(diag1->GetCellData()->GetArray("MatchingIdentifier"));
//    vtkIntArray* matchingIds2 = vtkIntArray::SafeDownCast(diag2->GetCellData()->GetArray("MatchingIdentifier"));
//    vtkIntArray* pairingIds1 = vtkIntArray::SafeDownCast(diag1->GetCellData()->GetArray("PairIdentifier"));
//    vtkIntArray* pairingIds2 = vtkIntArray::SafeDownCast(diag2->GetCellData()->GetArray("PairIdentifier"));
//    auto size1 = (int) diag1->GetCellData()->GetNumberOfTuples();
//    auto size2 = (int) diag2->GetCellData()->GetNumberOfTuples();
//
//    if (in % 2 == 0) {
//      for (int i = 0; i < size1; ++i) {
//        if (matchingIds1->GetValue(i) == -1) continue;
//        for (int j = 0; j < size2; ++j) {
//          if (matchingIds2->GetValue(j) == -1) continue;
//          if (matchingIds1->GetValue(i) == matchingIds2->GetValue(j))
//            matchings2.push_back(std::make_tuple(i, j, 0.0));
//        }
//      }
//    } else if (in % 2 == 1) {
//      for (int i = 0; i < size1; ++i) {
//        if (pairingIds1->GetValue(i) == -1) continue;
//        for (int j = 0; j < size2; ++j) {
//          if (pairingIds2->GetValue(j) == -1) continue;
//
//          if (pairingIds1->GetValue(i) == pairingIds2->GetValue(j))
//          {
//            matchings2.push_back(std::make_tuple(i, j, 0.0));
//          }
//        }
//      }
//    }
//
//    allMatchings.push_back(matchings2);
//    if (in == 0) continue;
//
//
//    // Do something
//    matchings1 = allMatchings[in - 1];
//
//    auto matchingsSize1 = (int) matchings1.size();
//    auto matchingsSize2 = (int) matchings2.size();
//    int endIndex = numPersistenceDiagramsInput - 2;
//
//    for (int i = 0; i < matchingsSize1; ++i) {
//      int m1ai0 = (int) std::get<0>(matchings1.at((unsigned long) i));
//      int m1ai1 = (int) std::get<1>(matchings1.at((unsigned long) i));
//
//      for (int j = 0; j < matchingsSize2; ++j) {
//        int m2aj0 = (int) std::get<0>(matchings2.at((unsigned long) j));
//        int m2aj1 = (int) std::get<1>(matchings2.at((unsigned long) j));
//
//        if (m1ai1 != m2aj0) continue;
//
//        // Detect in trackings and push.
//        bool found = false;
//        for (unsigned int k = 0; k < trackings.size(); ++k) {
//          trackingTuple *tt = &(trackings.at((unsigned long) k));
//          int chainStart = std::get<0>(*tt);
//          int chainEnd = std::get<1>(*tt);
//          std::vector<BIdVertex> chain = std::get<2>(*tt);
//
//          if (chainEnd == -1) {
//            auto chainSize = (int) chain.size();
//            if (chainSize == 0) {
//              // Should not happen
//              std::cout << "Brain error." << std::endl;
//            } else if (chainStart + chainSize == in &&
//              chain.at((unsigned long) chainSize - 1) == m1ai0) {
//              found = true;
//              chain.push_back(m1ai1);
//              int numEnd = in == endIndex ? endIndex : -1;
//              if (in == endIndex) {
//                chain.push_back(m2aj1);
//                std::get<1>(*tt) = numEnd;
//              }
//              std::get<2>(*tt) = chain;
//            }
//          }
//          (&trackings)->at((unsigned long) k) = std::make_tuple(chainStart, chainEnd, chain);
//        }
//
//        if (!found) {
//          std::vector<BIdVertex> chain;
//          chain.push_back(m1ai0);
//          chain.push_back(m1ai1);
//          if (in == endIndex) {
//            chain.push_back(m2aj1);
//          }
//          int numEnd = in == endIndex ? endIndex : -1;
//          trackingTuple tt = std::make_tuple(in - 1, numEnd, chain);
//          trackings.push_back(tt);
//        }
//        // Create new.
//      }
//    }
//
//    // End non-matched chains.
//    for (unsigned int k = 0; k < trackings.size(); ++k) {
//      trackingTuple *tt = &(trackings.at((unsigned long) k));
//      int chainStart = std::get<0>(*tt);
//      int chainEnd = std::get<1>(*tt);
//      if (chainEnd == -1) {
//        std::vector<BIdVertex> chain = std::get<2>(*tt);
//        int chainSize = (int) chain.size();
//        if (chainStart + chainSize - 1 < in)
//          std::get<1>(*tt) = in - 1;
//      }
//    }
//
//  }
//
//  std::sort(trackings.begin(), trackings.end(),
//    [](const trackingTuple &a, const trackingTuple &b) -> bool
//    {
//        return std::get<0>(a) < std::get<0>(b);
//    });
//
//  // Merge close connected components with threshold.
//  std::vector<std::set<int>> trackingTupleToMerged(trackings.size(), std::set<int>());
//  if (DoPostProc)
//    for (unsigned int k = 0; k < trackings.size(); ++k) {
//      trackingTuple tk = trackings.at((unsigned long) k);
//      int startK = std::get<0>(tk);
//      int endK = std::get<1>(tk);
//      if (endK < 0) endK = numPersistenceDiagramsInput - 1;
//      std::vector<BIdVertex> chainK = std::get<2>(tk);
//      std::vector<diagramTuple> diagramStartK = allDiagrams->at((unsigned long) startK);
//      std::vector<diagramTuple> diagramEndK = allDiagrams->at((unsigned long) endK);
//
//      int n1 = (int) chainK.at(0);
//      int n2 = (int) chainK.at(chainK.size() - 1);
//      diagramTuple tuple1 = diagramStartK.at((unsigned long) n1);
//      diagramTuple tuple2 = diagramEndK.at((unsigned long) n2);
//
//      double x1, y1, z1, x2, y2, z2;
//
//      BNodeType point1Type1 = std::get<1>(tuple1);
//      BNodeType point1Type2 = std::get<3>(tuple1);
//      bool t11Min = point1Type1 == BLocalMin; bool t11Max = point1Type1 == BLocalMax;
//      bool t12Min = point1Type2 == BLocalMin; bool t12Max = point1Type2 == BLocalMax;
//      // bool bothEx1 = t11Ex && t12Ex;
//      bool t1Max = t11Max || t12Max;
//      bool t1Min = !t1Max && (t11Min || t12Min);
//
//      x1 = t1Max ? std::get<11>(tuple1) : t1Min ? std::get<7>(tuple1) : 0;
//      y1 = t1Max ? std::get<12>(tuple1) : t1Min ? std::get<8>(tuple1) : 0;
//      z1 = t1Max ? std::get<13>(tuple1) : t1Min ? std::get<9>(tuple1) : 0;
//
//      BNodeType point2Type1 = std::get<1>(tuple2);
//      BNodeType point2Type2 = std::get<3>(tuple2);
//      bool t21Min = point2Type1 == BLocalMin; bool t21Max = point2Type1 == BLocalMax;
//      bool t22Min = point2Type2 == BLocalMin; bool t22Max = point2Type2 == BLocalMax;
//      // bool bothEx2 = t21Ex && t22Ex;
//      bool t2Max = t21Max || t22Max;
//      bool t2Min = !t2Max && (t21Min || t22Min);
//
//      // if (bothEx2) {
//      x2 = t2Max ? std::get<11>(tuple2) : t2Min ? std::get<7>(tuple2) : 0;
//      y2 = t2Max ? std::get<12>(tuple2) : t2Min ? std::get<8>(tuple2) : 0;
//      z2 = t2Max ? std::get<13>(tuple2) : t2Min ? std::get<9>(tuple2) : 0;
//      // }
//      // if (!bothEx1 && !bothEx2)
//      //  continue;
//
//      // Saddle-saddle matching not supported.
//      if (!t1Min && !t2Min && !t1Max && !t2Max)
//        continue;
//
//      // Check every other tracking trajectory.
//      for (unsigned int m = k + 1; m < trackings.size(); ++m) {
//        trackingTuple tm = trackings.at((unsigned long) m);
//        int startM = std::get<0>(tm);
//        int endM = std::get<1>(tm);
//        std::vector<BIdVertex> chainM = std::get<2>(tm);
//        if ((endK > 0 && startM > endK) || (endM > 0 && startK > endM)) continue;
//
//        for (int c = 0; c < (int) chainM.size(); ++c) {
//          bool doMatch1 = startM + c == startK;
//          bool doMatch2 = startM + c == endK;
//
//          // if (startM + c != startK && startM + c != endK) continue;
//          if (!doMatch1 && !doMatch2) continue;
//
//          /// Check proximity.
//          int n3 = (int) chainM.at(c);
//          std::vector<diagramTuple> diagramM = allDiagrams->at((unsigned long) startM + c);
//          diagramTuple tuple3 = diagramM.at((unsigned long) n3);
//          double x3, y3, z3;
//          BNodeType point3Type1 = std::get<1>(tuple3);
//          BNodeType point3Type2 = std::get<3>(tuple3);
//          bool t31Min = point3Type1 == BLocalMin; bool t31Max = point3Type1 == BLocalMax;
//          bool t32Min = point3Type2 == BLocalMin; bool t32Max = point3Type2 == BLocalMax;
//          // bool bothEx3 = t31Ex && t32Ex;
//          // if (!bothEx3)
//          //  continue;
//          bool t3Max = t31Max || t32Max;
//          bool t3Min = !t3Max && (t31Min || t32Min);
//
//          x3 = t3Max ? std::get<11>(tuple3) : t3Min ? std::get<7>(tuple3) : 0;
//          y3 = t3Max ? std::get<12>(tuple3) : t3Min ? std::get<8>(tuple3) : 0;
//          z3 = t3Max ? std::get<13>(tuple3) : t3Min ? std::get<9>(tuple3) : 0;
//
//          double dist = 0;
//          bool hasMatched = false;
//          if (doMatch1 && ((t3Max && t1Max) || (t3Min && t1Min))) {
//            double dist13 = sqrt(std::pow(x1 - x3, 2) + std::pow(y1 - y3, 2) + std::pow(z1 - z3, 2));
//            dist = dist13;
//            if (dist13 >= PostProcThresh) continue;
//            hasMatched = true;
//          }
//
//          if (doMatch2 && ((t3Max && t2Max) || (t3Min && t2Min))) {
//            double dist23 = sqrt(std::pow(x2 - x3, 2) + std::pow(y2 - y3, 2) + std::pow(z2 - z3, 2));
//            dist = dist23;
//            if (dist23 >= PostProcThresh) continue;
//            hasMatched = true;
//          }
//
//          if (!hasMatched)
//            continue;
//
//          /// Merge!
//          std::stringstream msg;
//          msg << "[ttkTrackingFromPersistenceDiagrams] Merged " << m << " with " << k << ": d = "
//              << dist << "." << std::endl;
//          dMsg(std::cout, msg.str(), timeMsg);
//
//          // Get every other tracking trajectory.
//          std::set<int>* mergedM = &(trackingTupleToMerged[m]);
//          // std::set<int> mergedK = trackingTupleToMerged[k];
//
//          // Push for others to merge.
//          // for (auto& i : mergedM) mergedK.insert(i);
//          // for (auto& i : mergedK) mergedM.insert(i);
//          // mergedK.insert(m);
//          mergedM->insert(k);
//          break;
//        }
//      }
//
//    }

  // Build mesh.
  int currentVertex = 0;
  for (unsigned int k = 0; k < trackings.size(); ++k) {
    trackingTuple tt = trackings.at((unsigned long) k);

    int numStart = std::get<0>(tt);
    // int numEnd = std::get<1>(tt);
    std::vector<BIdVertex> chain = std::get<2>(tt);

    if (chain.size() <= 1) {
      std::cout << "Got an unexpected 0-size chain." << std::endl;
      return -9;
    }

    for (int c = 0; c < (int) chain.size() - 1; c += 2)
    {
      if (numStart %2 != 0 && c == 0) c++;

      std::vector<matchingTuple> matchings1 = *(outputMatchings->at((unsigned long) numStart + c));
      int d1id = numStart + c;
      int d2id = d1id + 1; // c % 2 == 0 ? d1id + 1 : d1id;
      std::vector<diagramTuple> diagram1 = *(inputPersistenceDiagrams->at((unsigned long) d1id));
      std::vector<diagramTuple> diagram2 = *(inputPersistenceDiagrams->at((unsigned long) d2id));

      // Insert segments
      vtkIdType ids[2];
      int n1 = (int) chain.at((unsigned long) c);
      int n2 = (int) chain.at((unsigned long) c + 1);

      // Search for right matching.
      double cost = 0.0;
      for (int m = 0; m < (int) matchings1.size(); ++m) {
        matchingTuple tup = matchings1.at(m);
        int d1id1 = (int) std::get<0>(tup);
        if (d1id1 == n1) {
          cost = std::get<2>(tup);
          break;
        }
      }

      diagramTuple tuple1 = diagram1.at((unsigned long) n1);
      diagramTuple tuple2 = diagram2.at((unsigned long) n2);

      double x1, y1, z1, x2, y2, z2;

      BNodeType point1Type1 = std::get<1>(tuple1);
      BNodeType point1Type2 = std::get<3>(tuple1);
      BNodeType point1Type = point1Type1 == BLocalMax || point1Type2 == BLocalMax ? BLocalMax :
                             point1Type1 == BLocalMin || point1Type2 == BLocalMin ? BLocalMin :
                             point1Type1 == BSaddle2  || point1Type2 == BSaddle2  ? BSaddle2 : BSaddle1;
      bool t11Min = point1Type1 == BLocalMin; bool t11Max = point1Type1 == BLocalMax;
      bool t12Min = point1Type2 == BLocalMin; bool t12Max = point1Type2 == BLocalMax;
      bool bothEx1 = (t11Min && t12Max) || (t11Max && t12Min);
      if (bothEx1) {
        x1 = t12Max ? std::get<11>(tuple1) : std::get<7>(tuple1);
        y1 = t12Max ? std::get<12>(tuple1) : std::get<8>(tuple1);
        z1 = t12Max ? std::get<13>(tuple1) : std::get<9>(tuple1);
        if (useGeometricSpacing) z1 += spacing * (numStart + c / 2.0);
      } else {
        x1 = t12Max ? std::get<11>(tuple1) : t11Min ? std::get<7>(tuple1) : (std::get<7>(tuple1) + std::get<11>(tuple1)) / 2;
        y1 = t12Max ? std::get<12>(tuple1) : t11Min ? std::get<8>(tuple1) : (std::get<8>(tuple1) + std::get<12>(tuple1)) / 2;
        z1 = t12Max ? std::get<13>(tuple1) : t11Min ? std::get<9>(tuple1) : (std::get<9>(tuple1) + std::get<13>(tuple1)) / 2;
        if (useGeometricSpacing) z1 += spacing * (numStart + c / 2.0);
      }

      // Postproc component ids.
      int cid = k;
      bool hasMergedFirst = false;
      if (DoPostProc) {
        std::set<int> connected = trackingTupleToMerged->at(k);
        if (!connected.empty()) {
          int min = *(connected.begin());
          trackingTuple ttt = trackings.at((unsigned long) min);
          // int numStart2 = std::get<0>(ttt);
          int numEnd2 = std::get<1>(ttt);
          if ((numEnd2 > 0 && numStart + c > numEnd2 + 1) && min < (int) k) {
            // std::cout << "[ttkTrackingFromPersistenceDiagrams] Switched " << k << " for " << min << std::endl;
            cid = min;
            hasMergedFirst = numStart + c <= numEnd2 + 3;
          }

          if (hasMergedFirst) {
            std::cout << "Has merged first " << std::endl;

            // Replace former first end of the segment with previous ending segment.
            std::vector<BIdVertex> chain3 = std::get<2>(ttt);
            int nn = (int) chain3.at(chain3.size() - 1);
            std::vector<diagramTuple> diagramRematch = *(inputPersistenceDiagrams->at((unsigned long) numEnd2));
            diagramTuple tupleN = diagramRematch.at((unsigned long) nn);

            point1Type1 = std::get<1>(tupleN);
            point1Type2 = std::get<3>(tupleN);
            point1Type  = point1Type1 == BLocalMax || point1Type2 == BLocalMax ? BLocalMax :
                          point1Type1 == BLocalMin || point1Type2 == BLocalMin ? BLocalMin :
                          point1Type1 == BSaddle2  || point1Type2 == BSaddle2  ? BSaddle2 : BSaddle1;
            t11Min = point1Type1 == BLocalMin; t11Max = point1Type1 == BLocalMax;
            t12Min = point1Type2 == BLocalMin; t12Max = point1Type2 == BLocalMax;
            bothEx1 = (t11Min && t12Max) || (t11Max && t12Min);
            std::cout << "xyz " << x1 << ", " << y1 << ", " << z1 << std::endl;
            if (bothEx1) {
              x1 = t12Max ? std::get<11>(tupleN) : std::get<7>(tupleN);
              y1 = t12Max ? std::get<12>(tupleN) : std::get<8>(tupleN);
              z1 = t12Max ? std::get<13>(tupleN) : std::get<9>(tupleN);
              if (useGeometricSpacing) z1 += spacing * (numStart + c / 2.0);
            } else {
              x1 = t12Max ? std::get<11>(tupleN) : t11Min ? std::get<7>(tupleN) : (std::get<7>(tupleN) + std::get<11>(tupleN)) / 2;
              y1 = t12Max ? std::get<12>(tupleN) : t11Min ? std::get<8>(tupleN) : (std::get<8>(tupleN) + std::get<12>(tupleN)) / 2;
              z1 = t12Max ? std::get<13>(tupleN) : t11Min ? std::get<9>(tupleN) : (std::get<9>(tupleN) + std::get<13>(tupleN)) / 2;
              if (useGeometricSpacing) z1 += spacing * (numStart + c / 2.0);
            }
            std::cout << "xyz " << x1 << ", " << y1 << ", " << z1 << std::endl;
          }
        }
      }

      points->InsertNextPoint(x1, y1, z1);
      ids[0] = 2 * currentVertex;
      pointTypeScalars->InsertTuple1(ids[0], (double)(int) point1Type);
      timeScalars->InsertTuple1(ids[0], (double) numStart + c);
      componentIds->InsertTuple1(ids[0], cid);

      BNodeType point2Type1 = std::get<1>(tuple2);
      BNodeType point2Type2 = std::get<3>(tuple2);
      BNodeType point2Type = point2Type1 == BLocalMax || point2Type2 == BLocalMax ? BLocalMax :
                             point2Type1 == BLocalMin || point2Type2 == BLocalMin ? BLocalMin :
                             point2Type1 == BSaddle2 || point2Type2 == BSaddle2 ? BSaddle2 : BSaddle1;
      bool t21Ex = point2Type1 == BLocalMin || point2Type1 == BLocalMax;
      bool t22Ex = point2Type2 == BLocalMin || point2Type2 == BLocalMax;
      bool bothEx2 = t21Ex && t22Ex;
      if (bothEx2) {
        x2 = point2Type2 == BLocalMax ? std::get<11>(tuple2) : std::get<7>(tuple2);
        y2 = point2Type2 == BLocalMax ? std::get<12>(tuple2) : std::get<8>(tuple2);
        z2 = point2Type2 == BLocalMax ? std::get<13>(tuple2) : std::get<9>(tuple2);
        if (useGeometricSpacing) z2 += spacing * (numStart + c / 2.0 + 1);
      } else {
        x2 = t22Ex ? std::get<11>(tuple2) : t21Ex ? std::get<7>(tuple2) : (std::get<7>(tuple2) + std::get<11>(tuple2)) / 2;
        y2 = t22Ex ? std::get<12>(tuple2) : t21Ex ? std::get<8>(tuple2) : (std::get<8>(tuple2) + std::get<12>(tuple2)) / 2;
        z2 = t22Ex ? std::get<13>(tuple2) : t21Ex ? std::get<9>(tuple2) : (std::get<9>(tuple2) + std::get<13>(tuple2)) / 2;
        if (useGeometricSpacing) z2 += spacing * (numStart + c / 2.0 + 1);
      }
      points->InsertNextPoint(x2, y2, z2);
      ids[1] = 2 * currentVertex + 1;
      pointTypeScalars->InsertTuple1(ids[1], (double)(int) point2Type);
      timeScalars->InsertTuple1(ids[1], (double) numStart + c);

      // Postproc component ids.
      componentIds->InsertTuple1(ids[1], cid);

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);

      persistenceScalars->InsertTuple1(currentVertex, cost);
      valueScalars->InsertTuple1(currentVertex, (std::get<10>(tuple1) + std::get<10>(tuple2)) / 2);
      matchingIdScalars->InsertTuple1(currentVertex, currentVertex);

      currentVertex++;
    }

  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(valueScalars);
  persistenceDiagram->GetCellData()->AddArray(matchingIdScalars);
  persistenceDiagram->GetPointData()->AddArray(timeScalars);
  persistenceDiagram->GetPointData()->AddArray(componentIds);
  persistenceDiagram->GetPointData()->AddArray(pointTypeScalars);

  outputMesh_->ShallowCopy(persistenceDiagram);
  outputMesh->ShallowCopy(outputMesh_);

  std::cout << "Coucou, j'ai tout changé." << std::endl;
  delete outputPersistenceDiagrams;

  return 0;
}

int ttkTrackingFromPersistenceDiagrams::FillInputPortInformation(int port, vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkTrackingFromPersistenceDiagrams::FillOutputPortInformation(int port, vtkInformation *info) {
  if (port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}

int ttkTrackingFromPersistenceDiagrams::RequestData(vtkInformation *request,
  vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  ttk::Memory m;

  // Output pointers informations
  vtkInformation* outInfo;

  // Unified bound fields
  outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *mesh = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkUnstructuredGrid::DATA_OBJECT()));

  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  {
    std::stringstream msg;
    msg << "[ttkTrackingFromPersistenceDiagrams] Number of inputs: " << numInputs << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  // Get input data
  vtkDataSet* *input = new vtkDataSet*[numInputs];
  for (int i = 0; i < numInputs; i++)
  {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
  }

  doIt(input, mesh, numInputs);

  delete input;

  {
    std::stringstream msg;
    msg << "[ttkTrackingFromPersistenceDiagrams] Memory usage: " << m.getElapsedUsage()
      << " MB." << std::endl;
    dMsg(std::cout, msg.str(), memoryMsg);
  }

  return 1;
}
