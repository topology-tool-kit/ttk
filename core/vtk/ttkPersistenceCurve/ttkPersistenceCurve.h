/// \ingroup vtk
/// \class ttkPersistenceCurve
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date September 2016.
///
/// \brief TTK VTK-filter for the computation of persistence curves.
///
/// This filter computes the list of extremum-saddle pairs and computes the 
/// number of pairs as a function of persistence (i.e. the number of pairs 
/// whose persistence is higher than a threshold).
///
/// These curves provide useful visual clues in order to fine-tune persistence
/// simplification thresholds.
///
/// \param Input Input scalar field, either 2D or 3D, regular grid or 
/// triangulation (vtkDataSet)
/// \param Output0 Table giving the number of persistent minimum-saddle pairs 
/// as a function of persistence (vtkTable)
/// \param Output1 Table giving the number of persistent saddle-maximum pairs 
/// as a function of persistence (vtkTable)
/// \param Output2 Table giving the number of persistent all extremum-saddle 
/// pairs as a function of persistence (vtkTable)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttkFTMTreePP
/// \sa ttkPersistenceDiagram
/// \sa ttkTopologicalSimplification
/// \sa ttk::PersistenceCurve
#ifndef _TTK_PERSISTENCECURVE_H
#define _TTK_PERSISTENCECURVE_H

// ttk code includes
#include<PersistenceCurve.h>
#include<ttkWrapper.h>

// VTK includes
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkFiltersCoreModule.h>
#include<vtkInformation.h>
#include<vtkInformationVector.h>
#include<vtkDoubleArray.h>
#include<vtkTable.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>

class VTKFILTERSCORE_EXPORT ttkPersistenceCurve
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkPersistenceCurve* New();

    vtkTypeMacro(ttkPersistenceCurve, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }

    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);

    vtkSetMacro(UseInputOffsetScalarField, int);
    vtkGetMacro(UseInputOffsetScalarField, int);

    vtkSetMacro(InputOffsetScalarFieldName, string);
    vtkGetMacro(InputOffsetScalarFieldName, string);

    vtkSetMacro(ComputeSaddleConnectors, int);
    vtkGetMacro(ComputeSaddleConnectors, int);

    int getScalars(vtkDataSet* input);
    int getTriangulation(vtkDataSet* input);
    int getOffsets(vtkDataSet* input);

  protected:

    ttkPersistenceCurve();
    ~ttkPersistenceCurve();

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector, vtkInformationVector *outputVector);

    int FillOutputPortInformation(int port, vtkInformation* info);

    template <typename vtkArrayType, typename scalarType>
      int getPersistenceCurve(ftm::TreeType treeType, 
        const vector<pair<scalarType, ftm::idVertex>>& plot);

    template <typename vtkArrayType, typename scalarType>
      int getMSCPersistenceCurve(
        const vector<pair<scalarType, ftm::idVertex>>& plot);

  private:

    bool UseAllCores;
    int ThreadNumber;
    int ScalarFieldId;
    string ScalarField;
    string InputOffsetScalarFieldName;
    bool UseInputOffsetScalarField;
    bool ComputeSaddleConnectors;

    PersistenceCurve persistenceCurve_;
    Triangulation *triangulation_;
    vtkDataArray* inputScalars_;
    vtkTable* JTPersistenceCurve_;
    vtkTable* MSCPersistenceCurve_;
    vtkTable* STPersistenceCurve_;
    vtkTable* CTPersistenceCurve_;
    vtkIntArray* offsets_;
    vtkDataArray* inputOffsets_;
    bool varyingMesh_;
    vtkSmartPointer<ttkTriangulationFilter> inputTriangulation_;

    // base code features
    int doIt(vtkDataSet* input,
        vtkTable* outputJTPersistenceCurve,
        vtkTable* outputMSCPersistenceCurve,
        vtkTable* outputSTPersistenceCurve,
        vtkTable* outputCTPersistenceCurve);
    bool needsToAbort();
    int updateProgress(const float &progress);
};

template <typename vtkArrayType, typename scalarType>
int ttkPersistenceCurve::getPersistenceCurve(ftm::TreeType treeType, 
  const vector<pair<scalarType, ftm::idVertex>>& plot){
  const ftm::idVertex numberOfPairs = plot.size();

  vtkSmartPointer<vtkArrayType> persistenceScalars = 
    vtkSmartPointer<vtkArrayType>::New();
  vtkSmartPointer<vtkIntArray> numberOfPairsScalars = 
    vtkSmartPointer<vtkIntArray>::New();

  switch(treeType){
    case ftm::TreeType::Join:
      persistenceScalars->SetName("Persistence (minimum-saddle pairs)");
      numberOfPairsScalars->SetName("Number Of Pairs (minimum-saddle pairs)");
      break;

    case ftm::TreeType::Split:
      persistenceScalars->SetName("Persistence (maximum-saddle pairs)");
      numberOfPairsScalars->SetName("Number Of Pairs (maximum-saddle pairs)");
      break;

    case ftm::TreeType::Join_Split:
    case ftm::TreeType::Contour:
      persistenceScalars->SetName("Persistence (all pairs)");
      numberOfPairsScalars->SetName(
        "Number Of Pairs (all pairs)");
      break;
  }

  vtkSmartPointer<vtkTable> persistenceCurve = vtkSmartPointer<vtkTable>::New();

  if (numberOfPairs) {
    persistenceScalars->SetNumberOfTuples(numberOfPairs);
    numberOfPairsScalars->SetNumberOfTuples(numberOfPairs);
    persistenceCurve->SetNumberOfRows(numberOfPairs);

    for (ftm::idVertex i = 0; i < numberOfPairs; ++i) {
      persistenceScalars->SetTuple1(i, plot[i].first);
      numberOfPairsScalars->SetTuple1(i, plot[i].second);
    }

    persistenceCurve->AddColumn(persistenceScalars);
    persistenceCurve->AddColumn(numberOfPairsScalars);

    switch(treeType){
      case ftm::TreeType::Join:
        JTPersistenceCurve_->ShallowCopy(persistenceCurve);
        break;

      case ftm::TreeType::Split:
        STPersistenceCurve_->ShallowCopy(persistenceCurve);
        break;

      case ftm::TreeType::Join_Split:
      case ftm::TreeType::Contour:
        CTPersistenceCurve_->ShallowCopy(persistenceCurve);
        break;
    }
  }

  return 0;
}

template <typename vtkArrayType, typename scalarType>
int ttkPersistenceCurve::getMSCPersistenceCurve(
  const vector<pair<scalarType, ftm::idVertex>>& plot){
  const ftm::idVertex numberOfPairs = plot.size();

  vtkSmartPointer<vtkArrayType> persistenceScalars = 
    vtkSmartPointer<vtkArrayType>::New();
  vtkSmartPointer<vtkIntArray> numberOfPairsScalars = 
    vtkSmartPointer<vtkIntArray>::New();

  persistenceScalars->SetName("Persistence (saddle-saddle pairs)");
  numberOfPairsScalars->SetName("Number Of Pairs (saddle-saddle pairs)");

  vtkSmartPointer<vtkTable> persistenceCurve = vtkSmartPointer<vtkTable>::New();

  if (numberOfPairs) {
    persistenceScalars->SetNumberOfTuples(numberOfPairs);
    numberOfPairsScalars->SetNumberOfTuples(numberOfPairs);
    persistenceCurve->SetNumberOfRows(numberOfPairs);

    for (ftm::idVertex i = 0; i < numberOfPairs; ++i) {
      persistenceScalars->SetTuple1(i, plot[i].first);
      numberOfPairsScalars->SetTuple1(i, plot[i].second);
    }

    persistenceCurve->AddColumn(persistenceScalars);
    persistenceCurve->AddColumn(numberOfPairsScalars);

    MSCPersistenceCurve_->ShallowCopy(persistenceCurve);
  }

  return 0;
}

#endif // _TTK_PERSISTENCECURVE_H
