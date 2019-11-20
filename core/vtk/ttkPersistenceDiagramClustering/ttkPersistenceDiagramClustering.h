/// \ingroup base
/// \class ttk::ttkPersistenceDiagramBarycenter
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \brief TTK processing package for the computation of Wasserstein barycenters
/// and K-Means clusterings of a set of persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#ifndef _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H
#define _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H

#ifndef diagramTuple
#define diagramTuple                                                       \
  std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,            \
             ttk::CriticalType, dataType, ttk::SimplexId, dataType, float, \
             float, float, dataType, float, float, float>
#endif

#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// ttk code includes
#include <PersistenceDiagramClustering.h>
//
#include <PersistenceDiagramBarycenter.h>
//
#include <Wrapper.h>
//
#include <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceDiagramClustering
#else
class ttkPersistenceDiagramClustering
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  void setNumberOfInputsFromCommandLine(int number) {
    numberOfInputsFromCommandLine = number;
    SetNumberOfInputPorts(number);
  }
  static ttkPersistenceDiagramClustering *New();

  vtkTypeMacro(ttkPersistenceDiagramClustering, vtkDataSetAlgorithm);

  // default ttk setters
  void SetdebugLevel_(int debugLevel) {
    debugLevel_ = debugLevel;
    Modified();
    needUpdate_ = true;
  }

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
    needUpdate_ = true;
  }

  /*void SetThreadNumber(int threadNumber){
    ThreadNumber = threadNumber;
    SetThreads();
  }*/

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  // set-getters macros to define from each variable you want to access from
  // the outside (in particular from paraview) - to adapt.

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(WassersteinMetric, std::string);
  vtkGetMacro(WassersteinMetric, std::string);

  void SetUseProgressive(int data) {
    UseProgressive = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseProgressive, int);

  void SetTimeLimit(double data) {
    TimeLimit = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(TimeLimit, double);

  vtkSetMacro(UseOutputMatching, int);
  vtkGetMacro(UseOutputMatching, int);

  void SetThreadNumber(int data) {
    ThreadNumber = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(ThreadNumber, int);

  void SetAlpha(double data) {
    if(data > 0 && data <= 1) {

      Alpha = data;
    } else if(data > 1) {
      Alpha = 1;
    } else {
      Alpha = 0.001;
    }
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Alpha, double);

  void SetAntiAlpha(double data) {
    double alpha = 1 - data;
    SetAlpha(alpha);
  }

  void SetDeltaLim(double data) {
    DeltaLim = data;
    Modified();
    needUpdate_ = true;
  }

  vtkGetMacro(DeltaLim, double);

  void SetLambda(double data) {
    Lambda = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Lambda, double);

  void SetNumberOfClusters(int data) {
    NumberOfClusters = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(NumberOfClusters, int);

  void SetUseAccelerated(bool data) {
    UseAccelerated = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseAccelerated, bool);

  void SetUseKmeansppInit(bool data) {
    UseKmeansppInit = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseKmeansppInit, bool);

  void SetForceUseOfAlgorithm(bool data) {
    ForceUseOfAlgorithm = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(ForceUseOfAlgorithm, bool);

  void SetDeterministic(bool data) {
    Deterministic = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Deterministic, bool);

  void SetPairTypeClustering(int data) {
    PairTypeClustering = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(PairTypeClustering, int);

  void SetSpacing(double spacing) {
    Spacing = spacing;
    oldSpacing = spacing;
    Modified();
  }
  vtkGetMacro(Spacing, double);

  void SetDisplayMethod(int displayMethod) {
    DisplayMethod = displayMethod;
    if(displayMethod == 0) { // compact display
      Spacing = 0;
    } else {
      Spacing = oldSpacing;
    }
    Modified();
  }

  vtkGetMacro(DisplayMethod, bool);

  void SetUseAdditionalPrecision(bool data) {
    UseAdditionalPrecision = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseAdditionalPrecision, bool);

  void SetDistanceWritingOptions(int data) {
    DistanceWritingOptions = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(DistanceWritingOptions, int);

  void SetUseInterruptible(bool data) {
    UseInterruptible = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseInterruptible, bool);

  void SetMethod(int method) {
    Method = method;
    needUpdate_ = true;
    Modified();
  }
  vtkGetMacro(Method, double);

protected:
  ttkPersistenceDiagramClustering();

  ~ttkPersistenceDiagramClustering();

  template <typename dataType>
  double getPersistenceDiagram(std::vector<diagramTuple> *diagram,
                               vtkUnstructuredGrid *CTPersistenceDiagram_,
                               const double spacing,
                               const int diagramNumber);

  int FillInputPortInformation(int port, vtkInformation *info);
  int FillOutputPortInformation(int port, vtkInformation *info);

  template <typename dataType>
  vtkSmartPointer<vtkUnstructuredGrid>
    createMatchings(const vector<vector<diagramTuple>> *final_centroids,
                    vector<int> inv_clustering,
                    std::vector<std::vector<diagramTuple>> &all_CTDiagrams,
                    const vector<vector<vector<matchingTuple>>> *matchings,
                    double max_dimension,
                    double spacing);
  template <typename dataType>
  vtkSmartPointer<vtkUnstructuredGrid> createOutputClusteredDiagrams(
    std::vector<std::vector<diagramTuple>> &all_CTDiagrams,
    std::vector<int> inv_clustering,
    double max_dimension,
    double spacing);

  template <typename dataType>
  vtkSmartPointer<vtkUnstructuredGrid> createOutputCentroids(
    std::vector<std::vector<diagramTuple>> *final_centroids,
    std::vector<int> inv_clustering,
    double max_dimension,
    double spacing);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector);

private:
  // std::vector<std::vector<diagramTuple>>
  void *intermediateDiagrams_;
  // std::vector<std::vector<diagramTuple>>
  void *all_matchings_;
  void *final_centroids_;
  std::vector<int> inv_clustering_;

  // vtkUnstructuredGrid* output_clusters_;
  // vtkUnstructuredGrid* output_centroids_;

  int numberOfInputsFromCommandLine;
  int PairTypeClustering;
  bool ForceUseOfAlgorithm;
  bool Deterministic;
  bool UseAllCores;
  int ThreadNumber;
  bool UseOutputMatching;
  bool UseAdditionalPrecision;
  int DistanceWritingOptions;
  double Alpha;
  double DeltaLim;
  double Lambda;
  double Spacing;
  double oldSpacing;
  int DisplayMethod;
  bool UseInterruptible;
  int Method; // 0 = progressive approach, 1 = Auction approach
  double max_dimension_total_;

  bool needUpdate_;

  int NumberOfClusters;
  bool UseAccelerated;
  bool UseKmeansppInit;

  std::string ScalarField;
  std::string WassersteinMetric;

  bool UseProgressive;
  double TimeLimit;

  // base code features
  int doIt(vtkDataSet **input,
           vtkUnstructuredGrid *outputClusters,
           vtkUnstructuredGrid *outputCentroids,
           vtkUnstructuredGrid *outputMatchings,
           int numInputs);

  bool needsToAbort();

  int updateProgress(const float &progress);

  template <typename T>
  int dispatch(int numInputs,
               std::vector<vtkUnstructuredGrid *> &inputDiagram,
               vtkUnstructuredGrid *outputClusters,
               vtkUnstructuredGrid *outputCentroids,
               vtkUnstructuredGrid *outputMatchings);
};

template <typename dataType>
double ttkPersistenceDiagramClustering::getPersistenceDiagram(
  std::vector<diagramTuple> *diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber) {
  vtkIntArray *vertexIdentifierScalars
    = vtkIntArray::SafeDownCast(CTPersistenceDiagram_->GetPointData()->GetArray(
      ttk::VertexScalarFieldName));

  vtkIntArray *nodeTypeScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("CriticalType"));

  vtkIntArray *pairIdentifierScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray *extremumIndexScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairType"));

  vtkDoubleArray *persistenceScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("Persistence"));

  vtkDoubleArray *birthScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Birth"));
  vtkFloatArray *critCoordinates = vtkFloatArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Coordinates"));
  vtkDoubleArray *deathScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Death"));

  vtkPoints *points = (CTPersistenceDiagram_->GetPoints());
  int pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();
  // FIX : no more missed pairs
  for(int pair_index = 0; pair_index < pairingsSize; pair_index++) {
    const float index_of_pair = pair_index;
    if(*pairIdentifierScalars->GetTuple(pair_index) != -1)
      pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
  }
  // auto s = (float) 0.0;

  if(!deathScalars != !birthScalars)
    return -2;
  // bool Is3D = !(!deathScalars && !birthScalars);
  // if (!Is3D && diagramNumber == 1) s = (float) spacing;

  if(pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars
     || !nodeTypeScalars || !persistenceScalars || !extremumIndexScalars
     || !points)
    return -2;

  if(NumberOfClusters == 1) {
    diagram->resize(pairingsSize);
  } else {
    diagram->resize(pairingsSize + 1);
  }
  int nbNonCompact = 0;
  double max_dimension = 0;

  for(int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    double *critCoords1 = critCoordinates->GetTuple3(2 * i);

    auto coordX1 = (float)critCoords1[0];
    auto coordY1 = (float)critCoords1[1];
    auto coordZ1 = (float)critCoords1[2];

    double *critCoords2 = critCoordinates->GetTuple3(2 * i + 1);
    auto coordX2 = (float)critCoords2[0];
    auto coordY2 = (float)critCoords2[1];
    auto coordZ2 = (float)critCoords2[2];
    int index1 = 2 * i;
    double *coords1 = points->GetPoint(index1);
    auto x1 = (float)coords1[0];
    // auto y1 = (float) coords1[1];
    // auto z1 = (float) coords1[2];

    int index2 = index1 + 1;
    double *coords2 = points->GetPoint(index2);
    // auto x2 = (float) coords2[0];
    auto y2 = (float)coords2[1];
    // auto z2 = (float) coords2[2];

    dataType value1 = (!birthScalars) ? (dataType)x1
                                      : (dataType)birthScalars->GetValue(2 * i);
    dataType value2 = (!deathScalars)
                        ? (dataType)y2
                        : (dataType)deathScalars->GetValue(2 * i + 1);

    // if(value1 > max_dimension)  max_dimension = value1;
    // if(value2 > max_dimension)  max_dimension = value2;

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize) {
      if(pairIdentifier == 0) {
        max_dimension = (dataType)persistence;

        if(NumberOfClusters == 1) {
          diagram->at(0) = std::make_tuple(
            vertexId1, (BNodeType)0, vertexId2, (BNodeType)3,
            (dataType)persistence, pairType, value1, coordX1, coordY1, coordZ1,
            value2, coordX2, coordY2, coordZ2);
        } else {
          diagram->at(0) = std::make_tuple(
            vertexId1, (BNodeType)0, vertexId2, (BNodeType)1,
            (dataType)persistence, pairType, value1, coordX1, coordY1, coordZ1,
            value2, coordX2, coordY2, coordZ2);
          diagram->at(pairingsSize) = std::make_tuple(
            vertexId1, (BNodeType)1, vertexId2, (BNodeType)3,
            (dataType)persistence, pairType, value1, coordX1, coordY1, coordZ1,
            value2, coordX2, coordY2, coordZ2);
        }

      } else {
        diagram->at(pairIdentifier) = std::make_tuple(
          vertexId1, (BNodeType)nodeType1, vertexId2, (BNodeType)nodeType2,
          (dataType)persistence, pairType, value1, coordX1, coordY1, coordZ1,
          value2, coordX2, coordY2, coordZ2);
      }
    }
    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        std::stringstream msg;
        msg << "[TTKPersistenceDiagramClustering] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if(nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKPersistenceDiagramClustering] Missed " << nbNonCompact
          << " pairs due to non-compactness." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  return max_dimension;
}

template <typename dataType>
vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramClustering::createOutputCentroids(
    std::vector<std::vector<diagramTuple>> *final_centroids,
    std::vector<int> inv_clustering,
    double max_dimension,
    double spacing) {
  if(debugLevel_ > 5) {
    std::cout << "[ttkPersistenceDiagramClustering] Creating vtk diagrams"
              << std::endl;
  }
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkIntArray> nodeType = vtkSmartPointer<vtkIntArray>::New();
  nodeType->SetName("CriticalType");

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> idOfPair = vtkSmartPointer<vtkIntArray>::New();
  idOfPair->SetName("PairID");

  vtkSmartPointer<vtkDoubleArray> persistenceScalarsPoint
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalarsPoint->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> idOfDiagramPoint
    = vtkSmartPointer<vtkIntArray>::New();
  idOfDiagramPoint->SetName("ClusterID");

  vtkSmartPointer<vtkIntArray> pairType = vtkSmartPointer<vtkIntArray>::New();
  pairType->SetName("PairType");

  vtkSmartPointer<vtkFloatArray> coordsScalars
    = vtkSmartPointer<vtkFloatArray>::New();
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  int count = 0;
  for(unsigned int j = 0; j < final_centroids->size(); ++j) {
    std::vector<diagramTuple> *diagram = &((*final_centroids)[j]);

    // First, add diagram points to the global input diagram
    for(unsigned int i = 0; i < diagram->size(); ++i) {
      vtkIdType ids[2];
      diagramTuple t = diagram->at(i);
      double x1 = std::get<6>(t);
      double y1 = x1;
      if(DisplayMethod == 1 && spacing != 0) {
        x1 += 3 * (abs(spacing) + 0.2) * max_dimension * j;
      }
      double z1 = 0; // Change 1 to j if you want to isolate the diagrams

      float coords1[3];
      coords1[0] = std::get<7>(t);
      coords1[1] = std::get<8>(t);
      coords1[2] = std::get<9>(t);

      double x2 = std::get<6>(t);
      double y2 = std::get<10>(t);
      double z2 = 0; // Change 1 to j if you want to isolate the

      float coords2[3];
      coords2[0] = std::get<11>(t);
      coords2[1] = std::get<12>(t);
      coords2[2] = std::get<13>(t);

      idOfPair->InsertTuple1(count, i);

      points->InsertNextPoint(x1, y1, z1);
      coordsScalars->InsertTuple3(
        2 * count, coords1[0], coords1[1], coords1[2]);
      idOfDiagramPoint->InsertTuple1(2 * count, j);
      const ttk::CriticalType n1Type = std::get<1>(t);
      switch(n1Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count, 0);
      }
      if(DisplayMethod == 1 && spacing != 0) {
        points->InsertNextPoint(
          x2 + 3 * (abs(spacing) + 0.2) * max_dimension * j, y2, z2);
      } else {
        points->InsertNextPoint(x2, y2, z2);
      }
      coordsScalars->InsertTuple3(
        2 * count + 1, coords2[0], coords2[1], coords2[2]);
      idOfDiagramPoint->InsertTuple1(2 * count + 1, j);
      const ttk::CriticalType n2Type = std::get<3>(t);
      switch(n2Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count + 1, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count + 1, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count + 1, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count + 1, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count + 1, 0);
      }

      ids[0] = 2 * count;
      ids[1] = 2 * count + 1;

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      persistenceScalars->InsertTuple1(count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count + 1, y2 - x2);
      const ttk::SimplexId type = std::get<5>(t);
      switch(type) {
        case 0:
          pairType->InsertTuple1(count, 0);
          break;

        case 1:
          pairType->InsertTuple1(count, 1);
          break;

        case 2:
          pairType->InsertTuple1(count, 2);
          break;
        default:
          pairType->InsertTuple1(count, 0);
      }
      count++;
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(pairType);
  persistenceDiagram->GetCellData()->AddArray(idOfPair);
  persistenceDiagram->GetPointData()->AddArray(nodeType);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetPointData()->AddArray(idOfDiagramPoint);
  persistenceDiagram->GetPointData()->AddArray(persistenceScalarsPoint);

  return persistenceDiagram;
}

template <typename dataType>
vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramClustering::createOutputClusteredDiagrams(
    std::vector<std::vector<diagramTuple>> &all_CTDiagrams,
    std::vector<int> inv_clustering,
    double max_dimension,
    double spacing) {
  if(debugLevel_ > 5) {
    std::cout << "[ttkPersistenceDiagramClustering] Creating vtk Outputs"
              << std::endl;
  }
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkIntArray> nodeType = vtkSmartPointer<vtkIntArray>::New();
  nodeType->SetName("CriticalType");

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> idOfPair = vtkSmartPointer<vtkIntArray>::New();
  idOfPair->SetName("PairID");

  vtkSmartPointer<vtkDoubleArray> persistenceScalarsPoint
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalarsPoint->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> idOfDiagramPoint
    = vtkSmartPointer<vtkIntArray>::New();
  idOfDiagramPoint->SetName("DiagramID");

  vtkSmartPointer<vtkIntArray> idOfCluster
    = vtkSmartPointer<vtkIntArray>::New();
  idOfCluster->SetName("ClusterID");

  vtkSmartPointer<vtkIntArray> pairType = vtkSmartPointer<vtkIntArray>::New();
  pairType->SetName("PairType");

  vtkSmartPointer<vtkFloatArray> coordsScalars
    = vtkSmartPointer<vtkFloatArray>::New();
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  std::vector<int> cluster_size;
  std::vector<int> idxInCluster(all_CTDiagrams.size());
  for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j) {
    idxInCluster[j] = 0;
  }
  // RE-Invert clusters
  if(spacing > 0) {
    for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j) {
      unsigned int c = inv_clustering[j];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[j] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[j] = cluster_size[c] - 1;
      }
    }
  }
  int count = 0;
  for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j) {
    std::vector<diagramTuple> *diagram = &(all_CTDiagrams[j]);

    unsigned int c = inv_clustering[j];
    // First, add diagram points to the global input diagram
    for(unsigned int i = 0; i < diagram->size(); ++i) {
      vtkIdType ids[2];
      diagramTuple t = diagram->at(i);
      double x1 = std::get<6>(t);
      double y1 = x1;
      double z1 = 0;

      float coords1[3];
      coords1[0] = std::get<7>(t);
      coords1[1] = std::get<8>(t);
      coords1[2] = std::get<9>(t);
      double x2 = std::get<6>(t);
      double y2 = std::get<10>(t);
      double z2 = 0;
      if(DisplayMethod == 1 && spacing > 0) {
        // cout<<"j "<<j<<" size "<<cluster_size[inv_clustering[j]]<<endl;
        // cout<<"count "<<count_diagram<<endl;
        double angle = 2 * 3.1415926 * (double)(idxInCluster[j])
                       / cluster_size[inv_clustering[j]];
        x1 += (abs(spacing) + .2) * 3 * max_dimension * c
              + spacing * max_dimension * cos(angle);
        x2 += (abs(spacing) + .2) * 3 * max_dimension * c
              + spacing * max_dimension * cos(angle);
        y1 += spacing * max_dimension * sin(angle);
        y2 += spacing * max_dimension * sin(angle);
      } else if(DisplayMethod == 2) {
        z2 = spacing;
        z1 = spacing;
        if(j == 0) {
          z2 = -spacing;
          z1 = -spacing;
        }
      }

      float coords2[3];
      coords2[0] = std::get<11>(t);
      coords2[1] = std::get<12>(t);
      coords2[2] = std::get<13>(t);

      idOfPair->InsertTuple1(count, i);

      points->InsertNextPoint(x1, y1, z1);
      coordsScalars->InsertTuple3(
        2 * count, coords1[0], coords1[1], coords1[2]);
      idOfDiagramPoint->InsertTuple1(2 * count, j);
      // std::cout<<"\nMAX DIM \n"<<max_dimension<<std::endl;
      idOfCluster->InsertTuple1(2 * count, c);
      const ttk::CriticalType n1Type = std::get<1>(t);
      switch(n1Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count, 0);
      }

      points->InsertNextPoint(x2, y2, z2);
      coordsScalars->InsertTuple3(
        2 * count + 1, coords2[0], coords2[1], coords2[2]);
      idOfDiagramPoint->InsertTuple1(2 * count + 1, j);
      idOfCluster->InsertTuple1(2 * count + 1, c);
      const ttk::CriticalType n2Type = std::get<3>(t);
      switch(n2Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count + 1, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count + 1, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count + 1, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count + 1, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count + 1, 0);
      }

      ids[0] = 2 * count;
      ids[1] = 2 * count + 1;

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      persistenceScalars->InsertTuple1(count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count + 1, y2 - x2);
      const ttk::SimplexId type = std::get<5>(t);
      switch(type) {
        case 0:
          pairType->InsertTuple1(count, 0);
          break;

        case 1:
          pairType->InsertTuple1(count, 1);
          break;

        case 2:
          pairType->InsertTuple1(count, 2);
          break;
        default:
          pairType->InsertTuple1(count, 0);
      }
      count++;
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(pairType);
  persistenceDiagram->GetCellData()->AddArray(idOfPair);
  persistenceDiagram->GetPointData()->AddArray(nodeType);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetPointData()->AddArray(idOfDiagramPoint);
  persistenceDiagram->GetPointData()->AddArray(idOfCluster);
  persistenceDiagram->GetPointData()->AddArray(persistenceScalarsPoint);

  return persistenceDiagram;
}

template <typename dataType>
vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramClustering::createMatchings(
    const vector<vector<diagramTuple>> *final_centroids,
    vector<int> inv_clustering,
    std::vector<std::vector<diagramTuple>> &all_CTDiagrams,
    const vector<vector<vector<matchingTuple>>> *all_matchings,
    double max_dimension,
    double spacing) {
  if(debugLevel_ > 5) {
    std::cout << "[ttkPersistenceDiagramClustering] Creating vtk Matchings"
              << std::endl;
  }
  vtkSmartPointer<vtkPoints> matchingPoints = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnstructuredGrid> matchingMesh
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkIntArray> idOfDiagramMatchingPoint
    = vtkSmartPointer<vtkIntArray>::New();
  idOfDiagramMatchingPoint->SetName("DiagramID");

  vtkSmartPointer<vtkIntArray> idOfPoint = vtkSmartPointer<vtkIntArray>::New();
  idOfPoint->SetName("PointID");

  vtkSmartPointer<vtkIntArray> idOfDiagramMatching
    = vtkSmartPointer<vtkIntArray>::New();
  idOfDiagramMatching->SetName("DiagramID");

  vtkSmartPointer<vtkIntArray> idOfCluster
    = vtkSmartPointer<vtkIntArray>::New();
  idOfCluster->SetName("ClusterID");

  vtkSmartPointer<vtkDoubleArray> cost = vtkSmartPointer<vtkDoubleArray>::New();
  cost->SetName("Cost");

  vtkSmartPointer<vtkIntArray> pairType = vtkSmartPointer<vtkIntArray>::New();
  pairType->SetName("PairType");

  vtkSmartPointer<vtkIntArray> matchingCount
    = vtkSmartPointer<vtkIntArray>::New();
  matchingCount->SetName("MatchNumber");

  std::vector<int> cluster_size;
  std::vector<int> idxInCluster(all_CTDiagrams.size());

  std::vector<int> matchings_count(final_centroids->at(0).size(), 0);
  std::vector<int> count_to_good;

  for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j) {
    idxInCluster[j] = 0;
  }
  // RE-Invert clusters
  if(DisplayMethod == 1 && spacing > 0) {
    for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j) {
      unsigned int c = inv_clustering[j];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[j] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[j] = cluster_size[c] - 1;
      }
    }
  }
  int count = 0;
  for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j) {
    int c = inv_clustering[j];
    std::vector<diagramTuple> *diagram = &(all_CTDiagrams[j]);
    std::vector<matchingTuple> matchings_j
      = all_matchings->at(inv_clustering[j])[j];
    for(unsigned int i = 0; i < matchings_j.size(); ++i) {

      vtkIdType ids[2];
      ids[0] = 2 * count;
      ids[1] = 2 * count + 1;
      matchingTuple m = matchings_j[i];
      int bidder_id = std::get<0>(m);
      int good_id = std::get<1>(m);
      if(NumberOfClusters == 1 && good_id > -1) {
        matchings_count[good_id] += 1;
        count_to_good.push_back(good_id);
      }
      diagramTuple t1 = final_centroids->at(inv_clustering[j])[good_id];
      double x1 = std::get<6>(t1);
      double y1 = std::get<10>(t1);
      double z1 = 0;

      // endl;
      if(bidder_id < (int)diagram->size()) {
        diagramTuple t2 = diagram->at(bidder_id);
        double x2 = std::get<6>(t2);
        double y2 = std::get<10>(t2);
        double z2 = 0; // Change 1 to j if you want to isolate the diagrams

        if(DisplayMethod == 1 && spacing > 0) {
          double angle
            = 2 * 3.1415926 * (double)(idxInCluster[j]) / cluster_size[c];
          x1 += (abs(spacing) + .2) * 3 * max_dimension * c;
          x2 += (abs(spacing) + .2) * 3 * max_dimension * c
                + spacing * max_dimension * cos(angle);
          y2 += spacing * max_dimension * sin(angle);
        } else if(DisplayMethod == 2) {
          z2 = spacing;
          if(all_CTDiagrams.size() == 2 and j == 0) {
            z2 = -spacing;
          }
        }
        if(good_id > -1) {
          matchingPoints->InsertNextPoint(x1, y1, z1);
          matchingPoints->InsertNextPoint(x2, y2, z2);
          matchingMesh->InsertNextCell(VTK_LINE, 2, ids);
          idOfDiagramMatching->InsertTuple1(count, j);
          idOfCluster->InsertTuple1(count, inv_clustering[j]);
          cost->InsertTuple1(count, std::get<2>(m));
          idOfDiagramMatchingPoint->InsertTuple1(2 * count, j);
          idOfDiagramMatchingPoint->InsertTuple1(2 * count + 1, j);
          idOfPoint->InsertTuple1(2 * count, good_id);
          idOfPoint->InsertTuple1(2 * count + 1, bidder_id);

          const ttk::SimplexId type = std::get<5>(t2);
          switch(type) {
            case 0:
              pairType->InsertTuple1(count, 0);
              break;

            case 1:
              pairType->InsertTuple1(count, 1);
              break;

            case 2:
              pairType->InsertTuple1(count, 2);
              break;
            default:
              pairType->InsertTuple1(count, 0);
          }
          count++;
        }
      }
    }
  }
  if(NumberOfClusters == 1 and all_CTDiagrams.size() == 2) {
    for(int i = 0; i < count; i++) {
      matchingCount->InsertTuple1(i, matchings_count[count_to_good[i]]);
    }
  }

  matchingMesh->SetPoints(matchingPoints);
  matchingMesh->GetPointData()->AddArray(idOfDiagramMatchingPoint);
  matchingMesh->GetPointData()->AddArray(idOfPoint);
  matchingMesh->GetCellData()->AddArray(idOfDiagramMatching);
  matchingMesh->GetCellData()->AddArray(idOfCluster);
  matchingMesh->GetCellData()->AddArray(pairType);
  matchingMesh->GetCellData()->AddArray(cost);
  if(NumberOfClusters == 1 and all_CTDiagrams.size() == 2) {
    matchingMesh->GetCellData()->AddArray(matchingCount);
  }

  return matchingMesh;
}

#endif // _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H
