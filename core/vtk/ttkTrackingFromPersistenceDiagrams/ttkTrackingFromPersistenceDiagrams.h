#pragma once

#include <tuple>

#include <TrackingFromPersistenceDiagrams.h>
#include <ttkAlgorithm.h>
// #include <ttkUtils.h>
#include <ttkMacros.h>

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIndent.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkTrackingFromPersistenceDiagramsModule.h>

class TTKTRACKINGFROMPERSISTENCEDIAGRAMS_EXPORT
  ttkTrackingFromPersistenceDiagrams
  : public ttkAlgorithm,
    protected ttk::TrackingFromPersistenceDiagrams {

public:
  static ttkTrackingFromPersistenceDiagrams *New();

  vtkTypeMacro(ttkTrackingFromPersistenceDiagrams, ttkAlgorithm);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);

  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);

  vtkSetMacro(WassersteinMetric, const std::string &);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, const std::string &);
  vtkGetMacro(DistanceAlgorithm, std::string);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  vtkSetMacro(UseGeometricSpacing, bool);
  vtkGetMacro(UseGeometricSpacing, bool);

  vtkSetMacro(Spacing, double);
  vtkGetMacro(Spacing, double);

  vtkSetMacro(DoPostProc, bool);
  vtkGetMacro(DoPostProc, bool);

  vtkSetMacro(PostProcThresh, double);
  vtkGetMacro(PostProcThresh, double);

  using trackingTuple = ttk::trackingTuple;

  template <typename dataType>
  static int
    buildMesh(std::vector<trackingTuple> &trackings,
              std::vector<std::vector<ttk::MatchingType>> &outputMatchings,
              std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
              bool useGeometricSpacing,
              double spacing,
              bool DoPostProc,
              std::vector<std::set<int>> &trackingTupleToMerged,
              vtkPoints *points,
              vtkUnstructuredGrid *persistenceDiagram,
              vtkDoubleArray *persistenceScalars,
              vtkDoubleArray *valueScalars,
              vtkIntArray *matchingIdScalars,
              vtkIntArray *lengthScalars,
              vtkIntArray *timeScalars,
              vtkIntArray *componentIds,
              vtkIntArray *pointTypeScalars);

protected:
  ttkTrackingFromPersistenceDiagrams();

  ~ttkTrackingFromPersistenceDiagrams() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  // Warn: this is a duplicate from ttkBottleneckDistance.h
  template <typename dataType>
  int getPersistenceDiagram(
    ttk::DiagramType &diagram,
    const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram_,
    double spacing,
    int diagramNumber);

  // Warn: ditto
  template <typename dataType>
  int augmentPersistenceDiagrams(
    const ttk::DiagramType &diagram1,
    const ttk::DiagramType &diagram2,
    const std::vector<ttk::MatchingType> &matchings,
    const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram1_,
    const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram2_);

private:
  // Input bottleneck config.
  bool UseGeometricSpacing{false};
  bool Is3D{true};
  bool DoPostProc{false};
  double PostProcThresh{0.0};
  double Spacing{1.0};
  double Alpha{1.0};
  double Tolerance{1.0};
  double PX{1};
  double PY{1};
  double PZ{1};
  double PE{1};
  double PS{1};
  std::string DistanceAlgorithm{"ttk"};
  int PVAlgorithm{-1};
  std::string WassersteinMetric{"1"};

  vtkUnstructuredGrid *outputMesh_{nullptr};
};

template <typename dataType>
int ttkTrackingFromPersistenceDiagrams::buildMesh(
  std::vector<trackingTuple> &trackings,
  std::vector<std::vector<ttk::MatchingType>> &outputMatchings,
  std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
  bool useGeometricSpacing,
  double spacing,
  bool ttkNotUsed(DoPostProc),
  std::vector<std::set<int>> &trackingTupleToMerged,
  vtkPoints *points,
  vtkUnstructuredGrid *persistenceDiagram,
  vtkDoubleArray *persistenceScalars,
  vtkDoubleArray *valueScalars,
  vtkIntArray *matchingIdScalars,
  vtkIntArray *lengthScalars,
  vtkIntArray *timeScalars,
  vtkIntArray *componentIds,
  vtkIntArray *pointTypeScalars) {

  using ttk::CriticalType;
  int currentVertex = 0;
  for(unsigned int k = 0; k < trackings.size(); ++k) {
    trackingTuple tt = trackings.at((unsigned long)k);

    int numStart = std::get<0>(tt);
    //     int numEnd = std::get<1>(tt);
    std::vector<ttk::SimplexId> chain = std::get<2>(tt);
    int chainLength = chain.size();

    if(chain.size() <= 1) {
      // printErr("Got an unexpected 0-size chain.");
      return 0;
    }

    for(int c = 0; c < (int)chain.size() - 1; ++c) {
      std::vector<ttk::MatchingType> &matchings1
        = outputMatchings[numStart + c];
      int d1id = numStart + c;
      int d2id = d1id + 1; // c % 2 == 0 ? d1id + 1 : d1id;
      ttk::DiagramType &diagram1 = inputPersistenceDiagrams[d1id];
      ttk::DiagramType &diagram2 = inputPersistenceDiagrams[d2id];

      // Insert segments
      vtkIdType ids[2];
      auto n1 = (int)chain.at((unsigned long)c);
      auto n2 = (int)chain.at((unsigned long)c + 1);

      // Search for right matching.
      double cost = 0.0;
      for(const auto &tup : matchings1) {
        auto d1id1 = (int)std::get<0>(tup);
        if(d1id1 == n1) {
          cost = std::get<2>(tup);
          break;
        }
      }

      const auto &tuple1 = diagram1[n1];
      const auto &tuple2 = diagram2[n2];

      double x1, y1, z1, x2, y2, z2;

      auto point1Type1 = std::get<1>(tuple1);
      auto point1Type2 = std::get<3>(tuple1);
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
        x1 = t12Max ? std::get<11>(tuple1) : std::get<7>(tuple1);
        y1 = t12Max ? std::get<12>(tuple1) : std::get<8>(tuple1);
        z1 = t12Max ? std::get<13>(tuple1) : std::get<9>(tuple1);
        if(useGeometricSpacing)
          z1 += spacing * (numStart + c);
      } else {
        x1 = t12Max   ? std::get<11>(tuple1)
             : t11Min ? std::get<7>(tuple1)
                      : (std::get<7>(tuple1) + std::get<11>(tuple1)) / 2;
        y1 = t12Max   ? std::get<12>(tuple1)
             : t11Min ? std::get<8>(tuple1)
                      : (std::get<8>(tuple1) + std::get<12>(tuple1)) / 2;
        z1 = t12Max   ? std::get<13>(tuple1)
             : t11Min ? std::get<9>(tuple1)
                      : (std::get<9>(tuple1) + std::get<13>(tuple1)) / 2;
        if(useGeometricSpacing)
          z1 += spacing * (numStart + c);
      }

      // Postproc component ids.
      int cid = k;
      bool hasMergedFirst = false;
      // if(DoPostProc) {
      if(false) {
        std::set<int> &connected = trackingTupleToMerged[k];
        if(!connected.empty()) {
          int min = *(connected.begin());
          trackingTuple ttt = trackings.at((unsigned long)min);
          // int numStart2 = std::get<0>(ttt);
          int numEnd2 = std::get<1>(ttt);
          if((numEnd2 > 0 && numStart + c > numEnd2 + 1) && min < (int)k) {
            // std::cout << "[ttkTrackingFromPersistenceDiagrams] Switched " <<
            // k << " for " << min << std::endl;
            cid = min;
            hasMergedFirst = numStart + c <= numEnd2 + 3;
          }

          if(hasMergedFirst) {
            // std::cout << "Has merged first " << std::endl;

            // Replace former first end of the segment with previous ending
            // segment.
            std::vector<ttk::SimplexId> chain3 = std::get<2>(ttt);
            auto nn = (int)chain3.at(chain3.size() - 1);
            ttk::DiagramType &diagramRematch
              = inputPersistenceDiagrams[numEnd2];
            const auto &tupleN = diagramRematch.at((unsigned long)nn);

            point1Type1 = std::get<1>(tupleN);
            point1Type2 = std::get<3>(tupleN);
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
              x1 = t12Max ? std::get<11>(tupleN) : std::get<7>(tupleN);
              y1 = t12Max ? std::get<12>(tupleN) : std::get<8>(tupleN);
              z1 = t12Max ? std::get<13>(tupleN) : std::get<9>(tupleN);
              if(useGeometricSpacing)
                z1 += spacing * (numStart + c);
            } else {
              x1 = t12Max   ? std::get<11>(tupleN)
                   : t11Min ? std::get<7>(tupleN)
                            : (std::get<7>(tupleN) + std::get<11>(tupleN)) / 2;
              y1 = t12Max   ? std::get<12>(tupleN)
                   : t11Min ? std::get<8>(tupleN)
                            : (std::get<8>(tupleN) + std::get<12>(tupleN)) / 2;
              z1 = t12Max   ? std::get<13>(tupleN)
                   : t11Min ? std::get<9>(tupleN)
                            : (std::get<9>(tupleN) + std::get<13>(tupleN)) / 2;
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

      const auto point2Type1 = std::get<1>(tuple2);
      const auto point2Type2 = std::get<3>(tuple2);
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
      bool t21Ex = point2Type1 == CriticalType::Local_minimum
                   || point2Type1 == CriticalType::Local_maximum;
      bool t22Ex = point2Type2 == CriticalType::Local_minimum
                   || point2Type2 == CriticalType::Local_maximum;
      bool bothEx2 = t21Ex && t22Ex;
      if(bothEx2) {
        x2 = point2Type2 == CriticalType::Local_maximum ? std::get<11>(tuple2)
                                                        : std::get<7>(tuple2);
        y2 = point2Type2 == CriticalType::Local_maximum ? std::get<12>(tuple2)
                                                        : std::get<8>(tuple2);
        z2 = point2Type2 == CriticalType::Local_maximum ? std::get<13>(tuple2)
                                                        : std::get<9>(tuple2);
        if(useGeometricSpacing)
          z2 += spacing * (numStart + c + 1);
      } else {
        x2 = t22Ex   ? std::get<11>(tuple2)
             : t21Ex ? std::get<7>(tuple2)
                     : (std::get<7>(tuple2) + std::get<11>(tuple2)) / 2;
        y2 = t22Ex   ? std::get<12>(tuple2)
             : t21Ex ? std::get<8>(tuple2)
                     : (std::get<8>(tuple2) + std::get<12>(tuple2)) / 2;
        z2 = t22Ex   ? std::get<13>(tuple2)
             : t21Ex ? std::get<9>(tuple2)
                     : (std::get<9>(tuple2) + std::get<13>(tuple2)) / 2;
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
        currentVertex, (std::get<10>(tuple1) + std::get<10>(tuple2)) / 2.0);
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

// Warn: this is a duplicate from ttkBottleneckDistance
template <typename dataType>
int ttkTrackingFromPersistenceDiagrams::getPersistenceDiagram(
  ttk::DiagramType &diagram,
  const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber) {
  ttkSimplexIdTypeArray *vertexIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(
      CTPersistenceDiagram_->GetPointData()->GetArray(
        ttk::VertexScalarFieldName));

  vtkIntArray *nodeTypeScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("CriticalType"));

  ttkSimplexIdTypeArray *pairIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(
      CTPersistenceDiagram_->GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray *extremumIndexScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairType"));
  auto *persistenceScalars
    = CTPersistenceDiagram_->GetCellData()->GetArray("Persistence");
  auto *birthScalars = CTPersistenceDiagram_->GetCellData()->GetArray("Birth");

  const auto coords = vtkFloatArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray(
      ttk::PersistenceCoordinatesName));

  vtkPoints *points = (CTPersistenceDiagram_->GetPoints());
  if(!pairIdentifierScalars)
    return -2;

  auto pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();

  // Continuous indexing (no gap in indices)
  for(int pairIndex = 0; pairIndex < pairingsSize; ++pairIndex) {
    const float indexOfPair = pairIndex;
    if(*pairIdentifierScalars->GetTuple(pairIndex) != -1) // except diagonal
      pairIdentifierScalars->SetTuple(pairIndex, &indexOfPair);
  }

  auto s = (float)0.0;

  bool is2D = coords != nullptr;
  bool is3D = !is2D;
  if(Is3D && !is3D)
    Is3D = false;
  if(!Is3D && diagramNumber == 1)
    s = (float)spacing;

  if(pairingsSize < 1 || !vertexIdentifierScalars || !nodeTypeScalars
     || !persistenceScalars || !extremumIndexScalars || !points)
    return -2;

  diagram.resize((unsigned long)pairingsSize);
  int nbNonCompact = 0;

  for(int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetTuple1(i);

    int index1 = 2 * i;
    double *coords1 = points->GetPoint(index1);
    auto x1 = (float)coords1[0];
    auto y1 = (float)coords1[1];
    auto z1 = (float)coords1[2];

    int index2 = index1 + 1;
    double *coords2 = points->GetPoint(index2);
    auto x2 = (float)coords2[0];
    auto y2 = (float)coords2[1];
    auto z2 = (float)coords2[2];

    dataType value1 = birthScalars->GetTuple1(i);
    dataType value2 = value1 + persistence;

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize)
      diagram.at(pairIdentifier) = std::make_tuple(
        vertexId1, (ttk::CriticalType)nodeType1, vertexId2,
        (ttk::CriticalType)nodeType2, (dataType)persistence, pairType, value1,
        x1, y1, z1 + s, value2, x2, y2, z2 + s);

    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        this->printWrn("Diagram pair identifiers must be compact (not exceed "
                       "the diagram size).");
      }
    }
  }

  if(nbNonCompact > 0) {
    this->printWrn("Missed " + std::to_string(nbNonCompact)
                   + " pairs due to non-compactness.");
  }

  sort(diagram.begin(), diagram.end(),
       [](const ttk::PairTuple &a, const ttk::PairTuple &b) -> bool {
         return std::get<6>(a) < std::get<6>(b);
       });

  return 0;
}

// Warn: this is a duplicate from ttkBottleneckDistance
template <typename dataType>
int ttkTrackingFromPersistenceDiagrams::augmentPersistenceDiagrams(
  const ttk::DiagramType &diagram1,
  const ttk::DiagramType &diagram2,
  const std::vector<ttk::MatchingType> &matchings,
  const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram1_,
  const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram2_) {

  auto diagramSize1 = diagram1.size();
  auto diagramSize2 = diagram2.size();
  auto matchingsSize = matchings.size();

  vtkSmartPointer<vtkIntArray> matchingIdentifiers1
    = vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers1->SetName("MatchingIdentifier");

  vtkSmartPointer<vtkIntArray> matchingIdentifiers2
    = vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers2->SetName("MatchingIdentifier");

  if(matchingsSize > 0) {
    ttk::SimplexId ids[2];
    matchingIdentifiers1->SetNumberOfComponents(1);
    matchingIdentifiers2->SetNumberOfComponents(1);
    matchingIdentifiers1->SetNumberOfTuples(diagramSize1);
    matchingIdentifiers2->SetNumberOfTuples(diagramSize2);

    // Unaffected by default
    for(size_t i = 0; i < diagramSize1; ++i)
      matchingIdentifiers1->InsertTuple1(i, -1);
    for(size_t i = 0; i < diagramSize2; ++i)
      matchingIdentifiers2->InsertTuple1(i, -1);

    // Last cell = junction
    if(static_cast<int>(diagramSize1)
       < CTPersistenceDiagram1_->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers1->InsertTuple1(diagramSize1, -1);
      matchingIdentifiers1->InsertTuple1(diagramSize1 + 1, -1);
    }
    if(static_cast<int>(diagramSize2)
       < CTPersistenceDiagram2_->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers2->InsertTuple1(diagramSize2, -1);
      matchingIdentifiers2->InsertTuple1(diagramSize2 + 1, -1);
    }

    // Affect bottleneck matchings
    int pairingIndex = 0;
    for(size_t i = 0; i < matchingsSize; ++i) {
      const auto &t = matchings.at((unsigned long)i);
      ids[0] = std::get<0>(t);
      ids[1] = std::get<1>(t);
      matchingIdentifiers1->InsertTuple1(ids[0], pairingIndex);
      matchingIdentifiers2->InsertTuple1(ids[1], pairingIndex);
      pairingIndex++;
    }

    CTPersistenceDiagram1_->GetCellData()->AddArray(matchingIdentifiers1);
    CTPersistenceDiagram2_->GetCellData()->AddArray(matchingIdentifiers2);
  }

  return 0;
}
