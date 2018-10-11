/// \ingroup vtk
/// \class ttkPersistenceDiagramsClustering
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK VTK-filter that takes an input ensemble data set
/// (represented by a list of scalar fields) and which computes various
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \param Input0 Input ensemble scalar field #0 (vtkDataSet)
/// \param Input1 Input ensemble scalar field #1 (vtkDataSet)\n
/// ...\n
/// \param InputN Input ensemble scalar field #N (vtkDataSet)
/// \param Output0 Lower and upper bound fields (vtkDataSet)
/// \param Output1 Histogram estimations of the vertex probability density
/// functions (vtkDataSet)
/// \param Output2 Mean field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa vtkMandatoryCriticalPoints
/// \sa ttk::PersistenceDiagramsClustering
#ifndef _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H
#define _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H


#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1  ttk::CriticalType::Saddle1
#define BSaddle2  ttk::CriticalType::Saddle2

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkInformationVector.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkCellData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkTable.h>

// ttk code includes
#include                  <PersistenceDiagramsClustering.h>
#include                  <Wrapper.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceDiagramsClustering
#else
class ttkPersistenceDiagramsClustering
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkPersistenceDiagramsClustering* New();

    vtkTypeMacro(ttkPersistenceDiagramsClustering, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }

    /*void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }*/

    void SetUseAllCores(bool onOff){
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

	vtkSetMacro(UseProgressive, int);
    vtkGetMacro(UseProgressive, int);

    vtkSetMacro(TimeLimit, double);
    vtkGetMacro(TimeLimit, double);

	vtkSetMacro(UseOutputMatching, int);
    vtkGetMacro(UseOutputMatching, int);

	vtkSetMacro(ThreadNumber, int);
    vtkGetMacro(ThreadNumber, int);

	vtkSetMacro(Alpha, double);
    vtkGetMacro(Alpha, double);

	vtkSetMacro(NumberOfClusters, int);
	vtkGetMacro(NumberOfClusters, int);

	vtkSetMacro(UseAccelerated, bool);
	vtkGetMacro(UseAccelerated, bool);

	vtkSetMacro(UseKmeansppInit, bool);
	vtkGetMacro(UseKmeansppInit, bool);

  vtkSetMacro(Deterministic, bool);
  vtkGetMacro(Deterministic, bool);

  protected:

    ttkPersistenceDiagramsClustering();

    ~ttkPersistenceDiagramsClustering();

	template <typename dataType>
    int getPersistenceDiagram(
      std::vector<diagramTuple>* diagram,
      vtkUnstructuredGrid *CTPersistenceDiagram_,
      const double spacing,
      const int diagramNumber);

	int FillInputPortInformation(int port, vtkInformation *info);
    int FillOutputPortInformation(int port, vtkInformation *info);

    int RequestData(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);


  private:
    bool                  Deterministic;
    bool                  UseAllCores;
    int                   ThreadNumber;
	bool                  UseOutputMatching;
	double                Alpha;

	int 				  NumberOfClusters;
	bool 				  UseAccelerated;
	bool 				  UseKmeansppInit;

    std::string                ScalarField;
	std::string                WassersteinMetric;

	bool	              UseProgressive;
	double	              TimeLimit;

    // base code features
    int doIt(vtkDataSet **input,
             int numInputs);

    bool needsToAbort();

    int updateProgress(const float &progress);

};

template <typename dataType>
int ttkPersistenceDiagramsClustering::getPersistenceDiagram(
  std::vector<diagramTuple>* diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber)
{
  vtkIntArray* vertexIdentifierScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray(ttk::VertexScalarFieldName));

  vtkIntArray* nodeTypeScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("CriticalType"));

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
  bool Is3D = !(!deathScalars && !birthScalars);
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
      diagram->at(pairIdentifier) = std::make_tuple(
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
        std::stringstream msg;
        msg << "[TTKBottleneckDistance] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). "
            << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if (nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKBottleneckDistance] Missed "
          << nbNonCompact << " pairs due to non-compactness."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }


  return 0;
}

#endif // _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H
