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
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
/// example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///

#ifndef _TTK_PERSISTENCECURVE_H
#define _TTK_PERSISTENCECURVE_H
#pragma once

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
    vtkTable *outputCurve,
    ttk::ftm::TreeType treeType,
    const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot);

  template <typename vtkArrayType, typename scalarType>
  int getMSCPersistenceCurve(
    vtkTable *outputCurve,
    const std::vector<std::pair<scalarType, ttk::SimplexId>> &plot);

private:
  bool ForceInputOffsetScalarField{false};

  template <typename VTK_TT, typename TTK_TT>
  int dispatch(vtkTable *outputJTPersistenceCurve,
               vtkTable *outputMSCPersistenceCurve,
               vtkTable *outputSTPersistenceCurve,
               vtkTable *outputCTPersistenceCurve,
               const VTK_TT *inputScalars,
               const void *inputOffsets,
               const TTK_TT *triangulation);
};

#endif // _TTK_PERSISTENCECURVE_H
