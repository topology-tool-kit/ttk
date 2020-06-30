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
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

// VTK Module
#include <ttkPersistenceCurveModule.h>

// ttk code includes
#include <PersistenceCurve.h>
#include <ttkAlgorithm.h>

class TTKPERSISTENCECURVE_EXPORT ttkPersistenceCurve
  : public ttkAlgorithm,
    protected ttk::PersistenceCurve {

public:
  static ttkPersistenceCurve *New();

  vtkTypeMacro(ttkPersistenceCurve, ttkAlgorithm);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(ComputeSaddleConnectors, bool);
  vtkGetMacro(ComputeSaddleConnectors, bool);

  vtkTable *GetOutput();
  vtkTable *GetOutput(int);

protected:
  ttkPersistenceCurve();
  ~ttkPersistenceCurve() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  template <typename vtkArrayType, typename scalarType>
  int getPersistenceCurve(
    ttk::ftm::TreeType treeType,
    const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot);

  template <typename vtkArrayType, typename scalarType>
  int getMSCPersistenceCurve(
    const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot);

private:
  bool ForceInputOffsetScalarField{false};
  bool ComputeSaddleConnectors{false};

  template <typename VTK_TT, typename TTK_TT>
  int dispatch(std::vector<std::pair<scalarType, SimplexId>> &JTPlot,
    std::vector<std::pair<scalarType, SimplexId>> &STPlot,
    std::vector<std::pair<scalarType, SimplexId>> &MSCPlot,
    std::vector<std::pair<scalarType, SimplexId>> &CTPlot,
    const VTK_TT *inputScalars, const void *inputOffsets,
    const TTK_TT *triangulation);
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
