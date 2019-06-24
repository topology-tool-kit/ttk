/// \ingroup vtk
/// \class ttkBottleneckDistance
/// \author Maxime Soler <soler.maxime@total.com>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the bottleneckDistance processing package.
///
/// VTK wrapping code for the @BottleneckDistance package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::BottleneckDistance
#pragma once

// TTK includes
#include <BottleneckDistance.h>
#include <ttkWrapper.h>

// VTK includes
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// Misc.
#include <cstdlib>
#include <ctime>
#include <random>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkBottleneckDistance
#else
class ttkBottleneckDistance
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkBottleneckDistance *New();

  vtkTypeMacro(ttkBottleneckDistance, vtkDataSetAlgorithm);

  // Default ttk setters
  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

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

  vtkSetMacro(UseOutputMatching, int);
  vtkGetMacro(UseOutputMatching, int);

  vtkSetMacro(BenchmarkSize, int);
  vtkGetMacro(BenchmarkSize, int);

  vtkSetMacro(UsePersistenceMetric, int);
  vtkGetMacro(UsePersistenceMetric, int);

  vtkSetMacro(WassersteinMetric, std::string);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, std::string);
  vtkGetMacro(DistanceAlgorithm, std::string);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  vtkSetMacro(UseGeometricSpacing, int);
  vtkGetMacro(UseGeometricSpacing, int);

  vtkSetMacro(Spacing, double);
  vtkGetMacro(Spacing, double);

  vtkGetMacro(result, double);

  // Override input types.
  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      case 1:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      default:
        break;
    }
    return 1;
  }

  // Override output types.
  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      case 1:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      case 2:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      default:
        break;
    }
    return 1;
  }

  // Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
  template <typename dataType>
  int getPersistenceDiagram(
    std::vector<diagramTuple> &diagram,
    const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram_,
    double spacing,
    int diagramNumber);

  template <typename dataType>
  int generatePersistenceDiagram(std::vector<diagramTuple> &diagram, int size);

  // Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
  template <typename dataType>
  int augmentPersistenceDiagrams(
    const std::vector<diagramTuple> &diagram1,
    const std::vector<diagramTuple> &diagram2,
    const std::vector<matchingTuple> &matchings,
    vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram1_,
    vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram2_);

  template <typename dataType>
  int translateSecondDiagram(vtkUnstructuredGrid *&outputCT2, double &spacing);

  template <typename dataType>
  int getMatchingMesh(const std::vector<diagramTuple> &diagram1,
                      const std::vector<diagramTuple> &diagram2,
                      const std::vector<matchingTuple> &matchings,
                      bool useGeometricSpacing,
                      double spacing,
                      bool is2D);

  int doBenchmark();

protected:
  ttkBottleneckDistance() {
    // settings
    UseAllCores = false;
    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(3);

    // inputs
    BenchmarkSize = -1;
    UseOutputMatching = false;
    UsePersistenceMetric = false;
    UseGeometricSpacing = false;
    WassersteinMetric = "2";
    Alpha = 1.0;
    Tolerance = 1.0;
    PX = 0;
    PY = 0;
    PZ = 0;
    PE = 1;
    PS = 1;
    Spacing = 0.0;
    Is3D = false;
    PVAlgorithm = -1;

    // outputs
    result = -1.;
    CTPersistenceDiagram1_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    CTPersistenceDiagram2_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    CTPersistenceDiagram3_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
  }

  ~ttkBottleneckDistance(){};

  TTK_SETUP();

private:
  int BenchmarkSize;

  bool UseOutputMatching;
  bool Is3D;
  double Spacing;
  double Alpha;
  double Tolerance;
  double PX;
  double PY;
  double PZ;
  double PE;
  double PS;

  std::string DistanceAlgorithm;

  std::string WassersteinMetric;
  bool UsePersistenceMetric;
  bool UseGeometricSpacing;
  int PVAlgorithm;
  double result;

  vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram1_;
  vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram2_;
  vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram3_;

  ttk::BottleneckDistance bottleneckDistance_;
};

template <typename dataType>
int ttkBottleneckDistance::generatePersistenceDiagram(
  std::vector<diagramTuple> &diagram, const int size) {
  // srand(time(NULL));
  int vertexId1 = 1;
  int vertexId2 = 2;

  // diagram.resize(size);
  // diagram.push_back(std::make_tuple(
  //    0, BLocalMin,
  //    1, BLocalMax,
  //    1.0,
  //    2,
  //    0.0, 0.0001, 0.0001, 0.0,
  //    1.0, 0.0001, 1.0, 0.0
  //));

  std::default_random_engine generator1(1998985);
  std::default_random_engine generator2(8584584);
  std::uniform_real_distribution<> dis1(0.0, 1.0);
  std::uniform_real_distribution<> dis2(0.0, 1.0);

  for(int i = 1; i < size; ++i) {
    int r0 = 2; // (rand() % 3);
    float r1
      = 0.001f + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    // (float) dis1(generator1); // static_cast <float> (rand()) / static_cast
    // <float> (RAND_MAX);
    float r2
      = 0.001f + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    // (float) dis2(generator2); // static_cast <float> (rand()) / static_cast
    // <float> (RAND_MAX); r2 *= 0.1f;

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

    auto value1 = (dataType)x1;
    auto value2 = (dataType)y2;

    dataType persistence = y2 - x1;

    diagram.push_back(std::make_tuple(vertexId1, nodeType1, vertexId2,
                                      nodeType2, persistence, pairType, value1,
                                      x1, y1, z1, value2, x2, y2, z2));

    vertexId1++;
    vertexId2++;
  }

  sort(diagram.begin(), diagram.end(),
       [](const diagramTuple &a, const diagramTuple &b) -> bool {
         return std::get<6>(a) < std::get<6>(b);
       });

  return 0;
}

// Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
template <typename dataType>
int ttkBottleneckDistance::getPersistenceDiagram(
  std::vector<diagramTuple> &diagram,
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

  vtkDoubleArray *persistenceScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("Persistence"));

  vtkDoubleArray *birthScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Birth"));

  vtkDoubleArray *deathScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Death"));

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

  if(!deathScalars != !birthScalars)
    return -2;
  bool is2D = !deathScalars && !birthScalars;
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
    double persistence = persistenceScalars->GetValue(i);

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

    dataType value1 = (!birthScalars) ? (dataType)x1
                                      : (dataType)birthScalars->GetValue(2 * i);
    dataType value2 = (!deathScalars)
                        ? (dataType)y2
                        : (dataType)deathScalars->GetValue(2 * i + 1);

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize)
      diagram.at(pairIdentifier)
        = std::make_tuple(vertexId1, (BNodeType)nodeType1, vertexId2,
                          (BNodeType)nodeType2, (dataType)persistence, pairType,
                          value1, x1, y1, z1 + s, value2, x2, y2, z2 + s);

    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        std::stringstream msg;
        msg << "[TTKBottleneckDistance] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if(nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKBottleneckDistance] Missed " << nbNonCompact
          << " pairs due to non-compactness." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  sort(diagram.begin(), diagram.end(),
       [](const diagramTuple &a, const diagramTuple &b) -> bool {
         return std::get<6>(a) < std::get<6>(b);
       });

  return 0;
}

// Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
template <typename dataType>
int ttkBottleneckDistance::augmentPersistenceDiagrams(
  const std::vector<diagramTuple> &diagram1,
  const std::vector<diagramTuple> &diagram2,
  const std::vector<matchingTuple> &matchings,
  vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram1,
  vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram2) {

  auto diagramSize1 = (BIdVertex)diagram1.size();
  auto diagramSize2 = (BIdVertex)diagram2.size();
  auto matchingsSize = (BIdVertex)matchings.size();

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

  return 0;
}

template <typename dataType>
int ttkBottleneckDistance::translateSecondDiagram(
  vtkUnstructuredGrid *&outputCT2, double &spacing) {
  vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
  vtkPoints *points = (outputCT2->GetPoints());
  ttkSimplexIdTypeArray *pairIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(
      outputCT2->GetCellData()->GetArray("PairIdentifier"));
  auto pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();

  for(int i = 0; i < pairingsSize; ++i) {
    int index1 = 2 * i;
    double *coords1 = points->GetPoint(index1);
    auto x1 = coords1[0];
    auto y1 = coords1[1];
    auto z1 = coords1[2] + spacing;

    int index2 = index1 + 1;
    double *coords2 = points->GetPoint(index2);
    auto x2 = coords2[0];
    auto y2 = coords2[1];
    auto z2 = coords2[2] + spacing;
    points2->InsertNextPoint(x1, y1, z1);
    points2->InsertNextPoint(x2, y2, z2);
  }

  outputCT2->SetPoints(points2);

  return 0;
}

template <typename dataType>
int ttkBottleneckDistance::getMatchingMesh(
  const std::vector<diagramTuple> &diagram1,
  const std::vector<diagramTuple> &diagram2,
  const std::vector<matchingTuple> &matchings,
  const bool useGeometricSpacing,
  const double spacing,
  const bool is2D) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetName("Cost");

  vtkSmartPointer<vtkIntArray> matchingIdScalars
    = vtkSmartPointer<vtkIntArray>::New();
  matchingIdScalars->SetName("MatchingIdentifier");

  auto matchingsSize = (BIdVertex)matchings.size();

  // Build matchings.
  if(matchingsSize > 0) {

    for(BIdVertex i = 0; i < matchingsSize; ++i) {
      vtkIdType ids[2];

      matchingTuple t = matchings.at((unsigned long)i);
      auto n1 = (int)std::get<0>(t);
      auto n2 = (int)std::get<1>(t);

      diagramTuple tuple1 = diagram1.at(n1);
      diagramTuple tuple2 = diagram2.at(n2);

      double x1, y1, z1, x2, y2, z2;

      bool linkMiddles = false;
      if(linkMiddles) {
        x1 = (std::get<7>(tuple1) + std::get<11>(tuple1)) / 2;
        y1 = (std::get<8>(tuple1) + std::get<12>(tuple1)) / 2;
        z1 = (std::get<9>(tuple1) + std::get<13>(tuple1)) / 2;
      } else {
        BNodeType t11 = std::get<1>(tuple1);
        BNodeType t12 = std::get<3>(tuple1);
        bool t11Max = t11 == BLocalMin || t11 == BLocalMax;
        bool t12Max = t12 == BLocalMin || t12 == BLocalMax;
        if(is2D) { // Quickchage for highlighting 2D matching
          if(t11 != BLocalMax && t12 != BLocalMax) {
            t11Max = t11 != BLocalMin;
            t12Max = t12 != BLocalMin;
          }
        }
        x1 = t12Max ? std::get<11>(tuple1)
                    : t11Max ? std::get<7>(tuple1)
                             : (std::get<7>(tuple1) + std::get<11>(tuple1)) / 2;
        y1 = t12Max ? std::get<12>(tuple1)
                    : t11Max ? std::get<8>(tuple1)
                             : (std::get<8>(tuple1) + std::get<12>(tuple1)) / 2;
        z1 = t12Max ? std::get<13>(tuple1)
                    : t11Max ? std::get<9>(tuple1)
                             : (std::get<9>(tuple1) + std::get<13>(tuple1)) / 2;
      }
      points->InsertNextPoint(x1, y1, z1);

      if(linkMiddles) {
        x2 = (std::get<7>(tuple2) + std::get<11>(tuple2)) / 2;
        y2 = (std::get<8>(tuple2) + std::get<12>(tuple2)) / 2;
        z2 = (std::get<9>(tuple2) + std::get<13>(tuple2)) / 2;
        if(useGeometricSpacing)
          z2 += spacing;
      } else {
        BNodeType t21 = std::get<1>(tuple2);
        BNodeType t22 = std::get<3>(tuple2);
        bool t21Max = t21 == BLocalMin || t21 == BLocalMax;
        bool t22Max = t22 == BLocalMin || t22 == BLocalMax;
        if(is2D) { // Quickchage for highlighting 2D matching
          if(t21 != BLocalMax && t22 != BLocalMax) {
            t21Max = t21 != BLocalMin;
            t22Max = t22 != BLocalMin;
          }
        }
        x2 = t22Max ? std::get<11>(tuple2)
                    : t21Max ? std::get<7>(tuple2)
                             : (std::get<7>(tuple2) + std::get<11>(tuple2)) / 2;
        y2 = t22Max ? std::get<12>(tuple2)
                    : t21Max ? std::get<8>(tuple2)
                             : (std::get<8>(tuple2) + std::get<12>(tuple2)) / 2;
        z2 = t22Max ? std::get<13>(tuple2)
                    : t21Max ? std::get<9>(tuple2)
                             : (std::get<9>(tuple2) + std::get<13>(tuple2)) / 2;
      }
      points->InsertNextPoint(x2, y2, z2);

      ids[0] = 2 * i;
      ids[1] = 2 * i + 1;

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);

      persistenceScalars->InsertTuple1(i, std::get<2>(t));
      matchingIdScalars->InsertTuple1(i, i);
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(matchingIdScalars);

  CTPersistenceDiagram3_->ShallowCopy(persistenceDiagram);

  return 0;
}
