#include                  <ttkTrackingFromFields.h>

vtkStandardNewMacro(ttkTrackingFromFields)

constexpr unsigned int str2int(const char* str, int h = 0)
{
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

int ttkTrackingFromFields::FillOutputPortInformation(int port, vtkInformation* info)
{
  switch (port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    // default:
    //   break;
  }

  return 1;
}

int ttkTrackingFromFields::doIt(
  std::vector<vtkDataSet *> &inputs,
  std::vector<vtkDataSet *> &outputs)
{
  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  // triangulation_.setWrapper(this);
  // triangulation_.setInputData(input);
  internalTriangulation_ = ttkTriangulation::getTriangulation(input);
  internalTriangulation_->setWrapper(this);

  // Test validity of datasets (must present the same number of points).
  if (!input)
    return -1;

  // Get number and list of inputs.
  std::vector<vtkDataArray*> inputScalarFieldsRaw;
  std::vector<vtkDataArray*> inputScalarFields;
  int numberOfInputFields = input->GetPointData()->GetNumberOfArrays();
  if (numberOfInputFields < 3) {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] not enough input fields to perform tracking." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  vtkDataArray* firstScalarField = input->GetPointData()->GetArray(0);
  // inputScalarFields.push_back(firstScalarField);

  for (int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray* currentScalarField = input->GetPointData()->GetArray(i);
    if (!currentScalarField ||
        firstScalarField->GetDataType() != currentScalarField->GetDataType())
    {
      std::stringstream msg;
      msg << "[ttkTrackingFromField] inconsistent field data type or size (" << i << ")." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(
      inputScalarFieldsRaw.begin(),
      inputScalarFieldsRaw.end(),
      [] (vtkDataArray* a,
          vtkDataArray* b)
      {
          std::string s1 = a->GetName();
          std::string s2 = b->GetName();
          return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end());
      });

  int end = EndTimestep <= 0 ? numberOfInputFields : std::min(numberOfInputFields, EndTimestep);
  for (int i = StartTimestep; i < end; i += Sampling) {
    vtkDataArray* currentScalarField = inputScalarFieldsRaw[i];
    std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  // Input -> persistence filter.
  std::string algorithm = DistanceAlgorithm;
  int pvalg = PVAlgorithm;
  bool useTTKMethod = false;

  if (pvalg >= 0) {
    switch (pvalg) {
      case 0:
      case 1:
      case 2:
      case 3:
        useTTKMethod = true;
        break;
      case 4:
        break;
      default:
        std::stringstream msg;
        msg << "[ttkTrackingFromFieldsFromField] Unrecognized tracking method." << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  } else {
    switch (str2int(algorithm.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
      case str2int("1"):
      case str2int("legacy"):
      case str2int("2"):
      case str2int("geometric"):
      case str2int("3"):
      case str2int("parallel"):
        useTTKMethod = true;
        break;
      case str2int("4"):
      case str2int("greedy"):
        break;
      default:
        std::stringstream msg;
        msg << "[ttkTrackingFromField] Unrecognized tracking method." << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  }

  int res = 0;
  if (useTTKMethod)
    res = trackWithPersistenceMatching(
        input, output, inputScalarFields);
  else
  {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] The specified matching method does not." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return res;
}

// (*) Persistence-driven approach
int ttkTrackingFromFields::trackWithPersistenceMatching(
    vtkDataSet *input,
    vtkUnstructuredGrid *output,
    std::vector<vtkDataArray*> inputScalarFields)
{
  using dataType = double;
  unsigned long fieldNumber = inputScalarFields.size();

  // 0. get data
  trackingF_.setTriangulation(internalTriangulation_);
  std::vector<void *> inputFields(fieldNumber);
  for (int i = 0; i < (int) fieldNumber; ++i)
    inputFields[i] = inputScalarFields[i]->GetVoidPointer(0);
  trackingF_.setInputScalars(inputFields);

  // 0'. get offsets
  int numberOfVertices = (int) input->GetNumberOfPoints();
  vtkIdTypeArray* offsets_ = vtkIdTypeArray::New();
  offsets_->SetNumberOfComponents(1);
  offsets_->SetNumberOfTuples(numberOfVertices);
  offsets_->SetName("OffsetScalarField");
  for(int i = 0; i < numberOfVertices; ++i)
    offsets_->SetTuple1(i,i);
  trackingF_.setInputOffsets(offsets_->GetVoidPointer(0));

  // 1. get persistence diagrams.
  auto persistenceDiagrams = new std::vector<std::vector<diagramTuple>*>(fieldNumber);
//  std::vector<vtkUnstructuredGrid*> persistenceDiagrams(fieldNumber);

  trackingF_.performDiagramComputation((int) fieldNumber, persistenceDiagrams, this);

  // 2. call feature tracking with threshold.
  auto outputMatchings =
    new std::vector<std::vector<matchingTuple>*>(fieldNumber - 1);

  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = Is3D;
  std::string wasserstein = WassersteinMetric;

  tracking_.performMatchings(
      (int) fieldNumber,
      persistenceDiagrams,
      outputMatchings,
      algorithm, // Not from paraview, from enclosing tracking plugin
      wasserstein,
      tolerance,
      is3D,
      alpha, // Blending
      PX, PY, PZ, PS, PE, // Coefficients
      this // Wrapper for accessing threadNumber
  );

  for (int i = 0; i < (int) fieldNumber - 1; ++i)
  {
    vtkSmartPointer<ttkBottleneckDistance> bottleneckModule = vtkSmartPointer<ttkBottleneckDistance>::New();
  }

  outputMesh_ = vtkUnstructuredGrid::New();
  vtkUnstructuredGrid *outputMesh = vtkUnstructuredGrid::SafeDownCast(output);

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
  auto trackingsBase = new std::vector<trackingTuple>(); // trackings;
  tracking_.performTracking(
      persistenceDiagrams, outputMatchings,
      trackingsBase);

  auto trackingTupleToMerged = new std::vector<std::set<int>>(trackingsBase->size(), std::set<int>());
  if (DoPostProc)
    tracking_.performPostProcess(
      persistenceDiagrams,
      trackingsBase,
      trackingTupleToMerged,
      PostProcThresh);

  bool useGeometricSpacing = UseGeometricSpacing;
  std::vector<trackingTuple> trackings = *trackingsBase;

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
      std::vector<diagramTuple> diagram1 = *(persistenceDiagrams->at((unsigned long) d1id));
      std::vector<diagramTuple> diagram2 = *(persistenceDiagrams->at((unsigned long) d2id));

      // Insert segments
      vtkIdType ids[2];
      int n1 = (int) chain.at((unsigned long) c);
      int n2 = (int) chain.at((unsigned long) c + 1);

      // Search for right matching.
      double cost = 0.0;
      for (int m = 0; m < (int) matchings1.size(); ++m) {
        matchingTuple tup = matchings1.at((unsigned long) m);
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
            std::vector<diagramTuple> diagramRematch = *(persistenceDiagrams->at((unsigned long) numEnd2));
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

  std::cout << "Coucou, j'ai tout changé encore." << std::endl;

  return 0;
}
