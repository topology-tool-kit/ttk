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

// VTK includes
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

// ttk code includes
#include <PersistenceCurve.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceCurve
#else
class ttkPersistenceCurve
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkPersistenceCurve *New();

  vtkTypeMacro(ttkPersistenceCurve, vtkDataSetAlgorithm);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
  }

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ForceInputOffsetScalarField, int);
  vtkGetMacro(ForceInputOffsetScalarField, int);

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(ComputeSaddleConnectors, int);
  vtkGetMacro(ComputeSaddleConnectors, int);

  vtkTable *GetOutput();
  vtkTable *GetOutput(int);

  int getScalars(vtkDataSet *input);
  int getTriangulation(vtkDataSet *input);
  int getOffsets(vtkDataSet *input);

protected:
  ttkPersistenceCurve();
  ~ttkPersistenceCurve();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  template <typename vtkArrayType, typename scalarType>
  int getPersistenceCurve(
    ttk::ftm::TreeType treeType,
    const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot);

  template <typename vtkArrayType, typename scalarType>
  int getMSCPersistenceCurve(
    const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot);

private:
  bool UseAllCores;
  ttk::ThreadId ThreadNumber;
  int ScalarFieldId;
  int OffsetFieldId;
  std::string ScalarField;
  std::string InputOffsetScalarFieldName;
  bool ForceInputOffsetScalarField;
  bool ComputeSaddleConnectors;

  ttk::PersistenceCurve persistenceCurve_;
  ttk::Triangulation *triangulation_;
  vtkDataArray *inputScalars_;
  vtkTable *JTPersistenceCurve_;
  vtkTable *MSCPersistenceCurve_;
  vtkTable *STPersistenceCurve_;
  vtkTable *CTPersistenceCurve_;
  vtkDataArray *offsets_;
  vtkDataArray *inputOffsets_;
  bool varyingMesh_;
  vtkSmartPointer<ttkTriangulationFilter> inputTriangulation_;

  // base code features
  int doIt(vtkDataSet *input,
           vtkTable *outputJTPersistenceCurve,
           vtkTable *outputMSCPersistenceCurve,
           vtkTable *outputSTPersistenceCurve,
           vtkTable *outputCTPersistenceCurve);
  bool needsToAbort() override;
  int updateProgress(const float &progress) override;
};

template <typename vtkArrayType, typename scalarType>
int ttkPersistenceCurve::getPersistenceCurve(
  ttk::ftm::TreeType treeType,
  const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot) {
  const ttk::SimplexId numberOfPairs = plot.size();

  vtkSmartPointer<vtkArrayType> persistenceScalars
    = vtkSmartPointer<vtkArrayType>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> numberOfPairsScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  switch(treeType) {
    case ttk::ftm::TreeType::Join:
      persistenceScalars->SetName("Persistence (minimum-saddle pairs)");
      numberOfPairsScalars->SetName("Number Of Pairs (minimum-saddle pairs)");
      break;

    case ttk::ftm::TreeType::Split:
      persistenceScalars->SetName("Persistence (maximum-saddle pairs)");
      numberOfPairsScalars->SetName("Number Of Pairs (maximum-saddle pairs)");
      break;

    case ttk::ftm::TreeType::Join_Split:
    case ttk::ftm::TreeType::Contour:
      persistenceScalars->SetName("Persistence (all pairs)");
      numberOfPairsScalars->SetName("Number Of Pairs (all pairs)");
      break;
  }

  vtkSmartPointer<vtkTable> persistenceCurve = vtkSmartPointer<vtkTable>::New();

  if(numberOfPairs) {
    persistenceScalars->SetNumberOfTuples(numberOfPairs);
    numberOfPairsScalars->SetNumberOfTuples(numberOfPairs);
    persistenceCurve->SetNumberOfRows(numberOfPairs);

    for(ttk::SimplexId i = 0; i < numberOfPairs; ++i) {
      persistenceScalars->SetTuple1(i, plot[i].first);
      numberOfPairsScalars->SetTuple1(i, plot[i].second);
    }

    persistenceCurve->AddColumn(persistenceScalars);
    persistenceCurve->AddColumn(numberOfPairsScalars);

    switch(treeType) {
      case ttk::ftm::TreeType::Join:
        JTPersistenceCurve_->ShallowCopy(persistenceCurve);
        break;

      case ttk::ftm::TreeType::Split:
        STPersistenceCurve_->ShallowCopy(persistenceCurve);
        break;

      case ttk::ftm::TreeType::Join_Split:
      case ttk::ftm::TreeType::Contour:
        CTPersistenceCurve_->ShallowCopy(persistenceCurve);
        break;
    }
  }

  return 0;
}

template <typename vtkArrayType, typename scalarType>
int ttkPersistenceCurve::getMSCPersistenceCurve(
  const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot) {
  const ttk::SimplexId numberOfPairs = plot.size();

  vtkSmartPointer<vtkArrayType> persistenceScalars
    = vtkSmartPointer<vtkArrayType>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> numberOfPairsScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  persistenceScalars->SetName("Persistence (saddle-saddle pairs)");
  numberOfPairsScalars->SetName("Number Of Pairs (saddle-saddle pairs)");

  vtkSmartPointer<vtkTable> persistenceCurve = vtkSmartPointer<vtkTable>::New();

  if(numberOfPairs) {
    persistenceScalars->SetNumberOfTuples(numberOfPairs);
    numberOfPairsScalars->SetNumberOfTuples(numberOfPairs);
    persistenceCurve->SetNumberOfRows(numberOfPairs);

    for(ttk::SimplexId i = 0; i < numberOfPairs; ++i) {
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
