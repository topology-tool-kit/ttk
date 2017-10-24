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

// ttk code includes
#include                  <BottleneckDistance.h>
#include                  <ttkWrapper.h>

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkCellData.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
class VTKFILTERSCORE_EXPORT ttkBottleneckDistance 
  : public vtkDataSetAlgorithm, public Wrapper
{

  public:
      
    static ttkBottleneckDistance* New();
    
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

    vtkSetMacro(UseOutputMatching, int);
    vtkGetMacro(UseOutputMatching, int);

    vtkSetMacro(UsePersistenceMetric, int);
    vtkGetMacro(UsePersistenceMetric, int);

    vtkSetMacro(WassersteinMetric, string);
    vtkGetMacro(WassersteinMetric, string);

    vtkSetMacro(UseGeometricSpacing, int);
    vtkGetMacro(UseGeometricSpacing, int);

    vtkSetMacro(Spacing, double);
    vtkGetMacro(Spacing, double);

    vtkGetMacro(result, double);

    // Over-ride the input types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the input types, on a per input basis.
    // In this example, the first input type is forced to vtkUnstructuredGrid.
    // The second input type is forced to vtkImageData.
    int FillInputPortInformation(int port, vtkInformation *info) {

      switch (port) {
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

    // Over-ride the output types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the output types, on a per output basis.
    // In this example, the first output type is forced to vtkUnstructuredGrid.
    // The second output type is forced to vtkImageData.
    int FillOutputPortInformation(int port, vtkInformation *info) {

      switch (port) {
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

    template <typename scalarType>
    int getPersistenceDiagram(
      vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
        scalarType, float, float, float, scalarType, float, float, float>>* diagram,
      vtkUnstructuredGrid *CTPersistenceDiagram_, double spacing, int diagramNumber);

    template <typename scalarType>
    int augmentPersistenceDiagrams(
      vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
        scalarType, float, float, float, scalarType, float, float, float>>* diagram1,
      vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
        scalarType, float, float, float, scalarType, float, float, float>>* diagram2,
      vector<tuple<idVertex, idVertex, scalarType> >* matchings,
      vtkUnstructuredGrid *CTPersistenceDiagram1_,
      vtkUnstructuredGrid *CTPersistenceDiagram2_);

    template <typename scalarType>
    int getMatchingMesh(
      vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
        scalarType, float, float, float, scalarType, float, float, float>>* diagram1,
      vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
        scalarType, float, float, float, scalarType, float, float, float>>* diagram2,
      vector<tuple<idVertex, idVertex, scalarType> >* matchings,
      bool useGeometricSpacing, double spacing);

  protected:
   
    ttkBottleneckDistance() {
      
      // settings
      UseAllCores = false;
      SetNumberOfInputPorts(2);
      SetNumberOfOutputPorts(3);

      // inputs
      UsePersistenceMetric = false;
      UseGeometricSpacing = false;
      WassersteinMetric = "2";
      Alpha = 1.0;
      Spacing = 0.0;

      // outputs
      result = -1.;
      CTPersistenceDiagram3_ = vtkUnstructuredGrid::New();
    }
    
    ~ttkBottleneckDistance() {
      if(CTPersistenceDiagram3_)
        CTPersistenceDiagram3_->Delete();
    };
    
    TTK_SETUP();
    
  private:

    bool                  UseOutputMatching;
    bool                  Is3D;
    double                Spacing;
    double                Alpha;

    string                WassersteinMetric;
    bool                  UsePersistenceMetric;
    bool                  UseGeometricSpacing;
    double                result;

    vtkUnstructuredGrid*  CTPersistenceDiagram1_;
    vtkUnstructuredGrid*  CTPersistenceDiagram2_;
    vtkUnstructuredGrid*  CTPersistenceDiagram3_;

    BottleneckDistance    bottleneckDistance_;
    
};

template <typename scalarType>
int ttkBottleneckDistance::getPersistenceDiagram(
  vector<
      tuple<
          idVertex, NodeType, // vertex 1, type of vertex 1
          idVertex, NodeType, // vertex 2, type of vertex 2
          scalarType, idVertex, // persistence, coupling type
          scalarType, // function value at vertex 1
          float, float, float, // vertex 1 coordinates
          scalarType, // function value at vertex 2
          float, float, float // vertex 2 coordinates
      >
  >* diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  double spacing, int diagramNumber)
{

  vtkIntArray* vertexIdentifierScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->GetPointData()->GetArray("VertexIdentifier"));

  vtkIntArray* nodeTypeScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->GetPointData()->GetArray("NodeType"));

  vtkIntArray* pairIdentifierScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray* extremumIndexScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->GetCellData()->GetArray("PairType"));

  vtkDoubleArray* persistenceScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->GetCellData()->GetArray("Persistence"));

  vtkDoubleArray* birthScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->GetPointData()->GetArray("Birth"));

  vtkDoubleArray* deathScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->GetPointData()->GetArray("Death"));

  vtkPoints* points = (CTPersistenceDiagram_->GetPoints());
  int pairingsSize = (int) pairIdentifierScalars->GetNumberOfTuples();
  float s = (float) 0.0;
  if (!deathScalars && !birthScalars) {
    pairingsSize -= 2;
    Is3D = false;
    if (diagramNumber == 1)
      s = spacing;
  } else
    Is3D = true;

  if (pairingsSize < 1 || !vertexIdentifierScalars || !nodeTypeScalars ||
      !pairIdentifierScalars || !persistenceScalars || !extremumIndexScalars ||
      !points)
    return -2;

  diagram->resize(pairingsSize);
  for (int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2*i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2*i+1);
    int nodeType1 = nodeTypeScalars->GetValue(2*i);
    int nodeType2 = nodeTypeScalars->GetValue(2*i+1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    int index1 = 2*i;
    double* coords1 = points->GetPoint(index1);
    float x1 = (float) coords1[0];
    float y1 = (float) coords1[1];
    float z1 = (float) coords1[2];

    int index2 = index1 + 1;
    double* coords2 = points->GetPoint(index2);
    float x2 = (float) coords2[0];
    float y2 = (float) coords2[1];
    float z2 = (float) coords2[2];

    scalarType value1 = (!birthScalars) ? (scalarType) x1 : (scalarType) birthScalars->GetValue(2*i);
    scalarType value2 = (!deathScalars) ? (scalarType) y2 : (scalarType) deathScalars->GetValue(2*i+1);

    diagram->at(pairIdentifier) = make_tuple(
      vertexId1,
      (NodeType) nodeType1,
      vertexId2,
      (NodeType) nodeType2,
      (scalarType) persistence,
      pairType,
      value1,
      x1, y1, z1 + s,
      value2,
      x2, y2, z2 + s
    );
  }

  return 0;
}

template <typename scalarType>
int ttkBottleneckDistance::augmentPersistenceDiagrams(
  vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
    scalarType, float, float, float, scalarType, float, float, float>>* diagram1,
  vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
    scalarType, float, float, float, scalarType, float, float, float>>* diagram2,
  vector<tuple<idVertex, idVertex, scalarType> >* matchings,
  vtkUnstructuredGrid *CTPersistenceDiagram1_,
  vtkUnstructuredGrid *CTPersistenceDiagram2_)
{

  idVertex diagramSize1 = diagram1->size();
  idVertex diagramSize2 = diagram2->size();
  idVertex matchingsSize = matchings->size();

  vtkSmartPointer<vtkIntArray> matchingIdentifiers1 = vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers1->SetName("MatchingIdentifier");

  vtkSmartPointer<vtkIntArray> matchingIdentifiers2 = vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers2->SetName("MatchingIdentifier");

  if (matchingsSize > 0) {
    vtkIdType ids[2];
    matchingIdentifiers1->SetNumberOfComponents(1);
    matchingIdentifiers2->SetNumberOfComponents(1);

    // Unaffected by default
    for (idVertex i = 0; i < diagramSize1; ++i)
      matchingIdentifiers1->InsertTuple1(i, -1);
    for (idVertex i = 0; i < diagramSize2; ++i)
      matchingIdentifiers2->InsertTuple1(i, -1);
    // Last cell = junction
    if (diagramSize1 < CTPersistenceDiagram1_->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers1->InsertTuple1(diagramSize1, -1);
      matchingIdentifiers1->InsertTuple1(diagramSize1+1, -1);
    }
    if (diagramSize2 < CTPersistenceDiagram2_->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers2->InsertTuple1(diagramSize2, -1);
      matchingIdentifiers2->InsertTuple1(diagramSize2+1, -1);
    }

    // Affect bottleneck matchings
    int pairingIndex = 0;
    for (idVertex i = 0; i < matchingsSize; ++i) {
      tuple<idVertex, idVertex, scalarType> t = matchings->at(i);
      ids[0] = get<0>(t);
      ids[1] = get<1>(t);
      matchingIdentifiers1->InsertTuple1(ids[0], pairingIndex);
      matchingIdentifiers2->InsertTuple1(ids[1], pairingIndex);
      pairingIndex++;
    }

    CTPersistenceDiagram1_->GetCellData()->AddArray(matchingIdentifiers1);
    CTPersistenceDiagram2_->GetCellData()->AddArray(matchingIdentifiers2);
  }

  return 0;
}

template <typename scalarType>
int ttkBottleneckDistance::getMatchingMesh(
    vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
      scalarType, float, float, float, scalarType, float, float, float>>* diagram1,
    vector<tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
      scalarType, float, float, float, scalarType, float, float, float>>* diagram2,
    vector<tuple<idVertex, idVertex, scalarType> >* matchings,
    bool useGeometricSpacing, double spacing)
{
  
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> persistenceScalars =
      vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetName("Cost");

  vtkSmartPointer<vtkIntArray> matchingIdScalars =
      vtkSmartPointer<vtkIntArray>::New();
  matchingIdScalars->SetName("MatchingIdentifier");

  idVertex matchingsSize = matchings->size();

  // Build matchings.
  if (matchingsSize > 0) {

    for (idVertex i = 0; i < matchingsSize; ++i) {
      vtkIdType ids[2];

      tuple<idVertex, idVertex, scalarType> t = matchings->at(i);
      int n1 = get<0>(t);
      int n2 = get<1>(t);

      tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
        scalarType, float, float, float, scalarType, float, float, float> tuple1
        = diagram1->at(n1);

      tuple<idVertex, NodeType, idVertex, NodeType, scalarType, idVertex,
        scalarType, float, float, float, scalarType, float, float, float> tuple2
        = diagram2->at(n2);

      double x1, y1, z1, x2, y2, z2;

      if (Is3D) {
        x1 = (get<7>(tuple1) + get<11>(tuple1))/2;
        y1 = (get<8>(tuple1) + get<12>(tuple1))/2;
        z1 = (get<9>(tuple1) + get<13>(tuple1))/2;
      } else {
        NodeType t11 = get<1>(tuple1);
        NodeType t12 = get<3>(tuple1);
        bool t11Max = t11 == NodeType::Local_minimum || t11 == NodeType::Local_maximum;
        bool t12Max = t12 == NodeType::Local_minimum || t12 == NodeType::Local_maximum;
        x1 = t12Max ? get<11>(tuple1) :
             t11Max ? get<7>(tuple1) :
             (get<7>(tuple1) + get<11>(tuple1))/2;
        y1 = t12Max ? get<12>(tuple1) :
             t11Max ? get<8>(tuple1) :
             (get<8>(tuple1) + get<12>(tuple1))/2;
        z1 = t12Max ? get<13>(tuple1) :
             t11Max ? get<9>(tuple1) :
             (get<9>(tuple1) + get<13>(tuple1))/2;
      }
      points->InsertNextPoint(x1, y1, z1);

      if (Is3D) {
        x2 = (get<7>(tuple2) + get<11>(tuple2))/2;
        y2 = (get<8>(tuple2) + get<12>(tuple2))/2;
        z2 = (get<9>(tuple2) + get<13>(tuple2))/2;
        if (useGeometricSpacing) z2 += spacing;
      } else {
        NodeType t21 = get<1>(tuple2);
        NodeType t22 = get<3>(tuple2);
        bool t21Max = t21 == NodeType::Local_minimum || t21 == NodeType::Local_maximum;
        bool t22Max = t22 == NodeType::Local_minimum || t22 == NodeType::Local_maximum;

        x2 = t22Max ? get<11>(tuple2) :
             t21Max ? get<7>(tuple2) :
             (get<7>(tuple2) + get<11>(tuple2))/2;
        y2 = t22Max ? get<12>(tuple2) :
             t21Max ? get<8>(tuple2) :
             (get<8>(tuple2) + get<12>(tuple2))/2;
        z2 = t22Max ? get<13>(tuple2) :
             t21Max ? get<9>(tuple2) :
             (get<9>(tuple2) + get<13>(tuple2))/2;
      }
      points->InsertNextPoint(x2, y2, z2);

      ids[0] = 2*i;
      ids[1] = 2*i+1;

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);

      persistenceScalars->InsertTuple1(i, get<2>(t));
      matchingIdScalars->InsertTuple1(i, i);

    }

  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(matchingIdScalars);

  CTPersistenceDiagram3_->ShallowCopy(persistenceDiagram);

  return 0;
}
