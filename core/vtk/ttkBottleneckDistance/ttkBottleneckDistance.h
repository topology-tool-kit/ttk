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

// VTK includes
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

// TTK includes
#include                  <BottleneckDistance.h>
#include                  <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkBottleneckDistance
#else
class ttkBottleneckDistance
#endif
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

    // Override input types.
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

    // Override output types.
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

    template <typename dataType>
    int getPersistenceDiagram(
      vector<diagramTuple>* diagram,
      vtkUnstructuredGrid *CTPersistenceDiagram_,
      const double spacing,
      const int diagramNumber);

    template <typename dataType>
    int augmentPersistenceDiagrams(
      const vector<diagramTuple>* diagram1,
      const vector<diagramTuple>* diagram2,
      const vector<matchingTuple>* matchings,
      vtkUnstructuredGrid *CTPersistenceDiagram1_,
      vtkUnstructuredGrid *CTPersistenceDiagram2_);

    template <typename dataType>
    int getMatchingMesh(
      const vector<diagramTuple>* diagram1,
      const vector<diagramTuple>* diagram2,
      const vector<matchingTuple>* matchings,
      const bool useGeometricSpacing,
      const double spacing);

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

template <typename dataType>
int ttkBottleneckDistance::getPersistenceDiagram(
  vector<diagramTuple>* diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber)
{
  vtkIntArray* vertexIdentifierScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("VertexIdentifier"));

  vtkIntArray* nodeTypeScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("NodeType"));

  vtkIntArray* pairIdentifierScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray* extremumIndexScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetCellData()->GetArray("PairType"));

  vtkDoubleArray* persistenceScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
      GetCellData()->GetArray("Persistence"));

  vtkDoubleArray* birthScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("Birth"));

  vtkDoubleArray* deathScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("Death"));

  vtkPoints* points = (CTPersistenceDiagram_->GetPoints());
  auto pairingsSize = (int) pairIdentifierScalars->GetNumberOfTuples();
  auto s = (float) 0.0;

  if (!deathScalars != !birthScalars) return -2;
  Is3D = !(!deathScalars && !birthScalars);
  if (!Is3D && diagramNumber == 1) s = (float) spacing;

  if (pairingsSize < 1 || !vertexIdentifierScalars
      || !pairIdentifierScalars || !nodeTypeScalars
      || !persistenceScalars || !extremumIndexScalars || !points)
    return -2;

  diagram->resize(pairingsSize);
  int nbNonCompact = 0;

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
    auto x1 = (float) coords1[0];
    auto y1 = (float) coords1[1];
    auto z1 = (float) coords1[2];

    int index2 = index1 + 1;
    double* coords2 = points->GetPoint(index2);
    auto x2 = (float) coords2[0];
    auto y2 = (float) coords2[1];
    auto z2 = (float) coords2[2];

    dataType value1 = (!birthScalars) ? (dataType) x1 :
                      (dataType) birthScalars->GetValue(2*i);
    dataType value2 = (!deathScalars) ? (dataType) y2 :
                      (dataType) deathScalars->GetValue(2*i+1);

    if (pairIdentifier != -1 && pairIdentifier < pairingsSize)
      diagram->at(pairIdentifier) = make_tuple(
        vertexId1, (BNodeType) nodeType1,
        vertexId2, (BNodeType) nodeType2,
        (dataType) persistence,
        pairType,
        value1, x1, y1, z1 + s,
        value2, x2, y2, z2 + s
      );

    if (pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if (nbNonCompact == 0) {
        stringstream msg;
        msg << "[TTKBottleneckDistance] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). "
            << endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if (nbNonCompact > 0) {
    {
      stringstream msg;
      msg << "[TTKBottleneckDistance] Missed "
          << nbNonCompact << " pairs due to non-compactness."
          << endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }


  return 0;
}

template <typename dataType>
int ttkBottleneckDistance::augmentPersistenceDiagrams(
  const vector<diagramTuple>* diagram1,
  const vector<diagramTuple>* diagram2,
  const vector<matchingTuple>* matchings,
  vtkUnstructuredGrid *CTPersistenceDiagram1_,
  vtkUnstructuredGrid *CTPersistenceDiagram2_)
{

  BIdVertex diagramSize1 = diagram1->size();
  BIdVertex diagramSize2 = diagram2->size();
  BIdVertex matchingsSize = matchings->size();

  vtkSmartPointer<vtkIntArray> matchingIdentifiers1 =
    vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers1->SetName("MatchingIdentifier");

  vtkSmartPointer<vtkIntArray> matchingIdentifiers2 =
    vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers2->SetName("MatchingIdentifier");

  if (matchingsSize > 0) {
    vtkIdType ids[2];
    matchingIdentifiers1->SetNumberOfComponents(1);
    matchingIdentifiers2->SetNumberOfComponents(1);

    // Unaffected by default
    for (BIdVertex i = 0; i < diagramSize1; ++i)
      matchingIdentifiers1->InsertTuple1(i, -1);
    for (BIdVertex i = 0; i < diagramSize2; ++i)
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
    for (BIdVertex i = 0; i < matchingsSize; ++i) {
      matchingTuple t = matchings->at(i);
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

template <typename dataType>
int ttkBottleneckDistance::getMatchingMesh(
    const vector<diagramTuple>* diagram1,
    const vector<diagramTuple>* diagram2,
    const vector<matchingTuple>* matchings,
    const bool useGeometricSpacing,
    const double spacing)
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

  BIdVertex matchingsSize = matchings->size();

  // Build matchings.
  if (matchingsSize > 0) {

    for (BIdVertex i = 0; i < matchingsSize; ++i) {
      vtkIdType ids[2];

      matchingTuple t = matchings->at(i);
      int n1 = get<0>(t);
      int n2 = get<1>(t);

      diagramTuple tuple1 = diagram1->at(n1);
      diagramTuple tuple2 = diagram2->at(n2);

      double x1, y1, z1, x2, y2, z2;

      if (Is3D) {
        x1 = (get<7>(tuple1) + get<11>(tuple1))/2;
        y1 = (get<8>(tuple1) + get<12>(tuple1))/2;
        z1 = (get<9>(tuple1) + get<13>(tuple1))/2;
      } else {
        BNodeType t11 = get<1>(tuple1);
        BNodeType t12 = get<3>(tuple1);
        bool t11Max = t11 == BLocalMin || t11 == BLocalMax;
        bool t12Max = t12 == BLocalMin || t12 == BLocalMax;
        x1 = t12Max ? get<11>(tuple1) : t11Max ? get<7>(tuple1) :
             (get<7>(tuple1) + get<11>(tuple1))/2;
        y1 = t12Max ? get<12>(tuple1) : t11Max ? get<8>(tuple1) :
             (get<8>(tuple1) + get<12>(tuple1))/2;
        z1 = t12Max ? get<13>(tuple1) : t11Max ? get<9>(tuple1) :
             (get<9>(tuple1) + get<13>(tuple1))/2;
      }
      points->InsertNextPoint(x1, y1, z1);

      if (Is3D) {
        x2 = (get<7>(tuple2) + get<11>(tuple2))/2;
        y2 = (get<8>(tuple2) + get<12>(tuple2))/2;
        z2 = (get<9>(tuple2) + get<13>(tuple2))/2;
        if (useGeometricSpacing) z2 += spacing;
      } else {
        BNodeType t21 = get<1>(tuple2);
        BNodeType t22 = get<3>(tuple2);
        bool t21Max = t21 == BLocalMin || t21 == BLocalMax;
        bool t22Max = t22 == BLocalMin || t22 == BLocalMax;

        x2 = t22Max ? get<11>(tuple2) : t21Max ? get<7>(tuple2) :
             (get<7>(tuple2) + get<11>(tuple2))/2;
        y2 = t22Max ? get<12>(tuple2) : t21Max ? get<8>(tuple2) :
             (get<8>(tuple2) + get<12>(tuple2))/2;
        z2 = t22Max ? get<13>(tuple2) : t21Max ? get<9>(tuple2) :
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
