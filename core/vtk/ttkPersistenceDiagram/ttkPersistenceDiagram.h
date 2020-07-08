/// \ingroup vtk
/// \class ttkPersistenceDiagram
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2016.
///
/// \brief TTK VTK-filter for the computation of persistence diagrams.
///
/// This filter computes the persistence diagram of the extremum-saddle pairs
/// of an input scalar field. The X-coordinate of each pair corresponds to its
/// birth, while its smallest and highest Y-coordinates correspond to its birth
/// and death respectively.
///
/// In practice, the diagram is represented by a vtkUnstructuredGrid. Each
/// vertex of this mesh represent a critical point of the input data. It is
/// associated with point data (vertexId, critical type). Each vertical edge
/// of this mesh represent a persistence pair. It is associated with cell data
/// (persistence of the pair, critical index of the extremum of the pair).
///
/// Persistence diagrams are useful and stable concise representations of the
/// topological features of a data-set. It is useful to fine-tune persistence
/// thresholds for topological simplification or for fast similarity
/// estimations for instance.
///
/// \param Input Input scalar field, either 2D or 3D, regular grid or
/// triangulation (vtkDataSet)
/// \param Output Output persistence diagram (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Computational Topology: An Introduction" \n
/// Herbert Edelsbrunner and John Harer \n
/// American Mathematical Society, 2010
///
/// \sa ttkFTMTreePP
/// \sa ttkPersistenceCurve
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkTopologicalSimplification
/// \sa ttk::PersistenceDiagram
#ifndef _TTK_PERSISTENCEDIAGRAM_H
#define _TTK_PERSISTENCEDIAGRAM_H
#pragma once

// VTK includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceDiagramModule.h>

// ttk code includes
#include <PersistenceDiagram.h>
#include <ttkAlgorithm.h>

class TTKPERSISTENCEDIAGRAM_EXPORT ttkPersistenceDiagram
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagram {

public:
  static ttkPersistenceDiagram *New();

  vtkTypeMacro(ttkPersistenceDiagram, ttkAlgorithm);

  void SetForceInputOffsetScalarField(int data) {
    ForceInputOffsetScalarField = data;
    Modified();
    computeDiagram_ = true;
  }
  vtkGetMacro(ForceInputOffsetScalarField, int);

  void SetComputeSaddleConnectors(int data) {
    ComputeSaddleConnectors = data;
    Modified();
    computeDiagram_ = true;
  }
  vtkGetMacro(ComputeSaddleConnectors, int);

  void SetShowInsideDomain(int onOff) {
    ShowInsideDomain = onOff;
    Modified();
  }
  vtkGetMacro(ShowInsideDomain, int);

  template <typename scalarType, typename vtkSimplexArray, class triangulationType>
  int setPersistenceDiagramInfo(
    ttk::SimplexId id,
    vtkSmartPointer<vtkSimplexArray> vertexIdentifierScalars,
    vtkSmartPointer<vtkIntArray> nodeTypeScalars,
    vtkSmartPointer<vtkFloatArray> coordsScalars,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram,
    vtkSmartPointer<vtkPoints> points,
    vtkIdType ids[3],
    vtkDataArray *inputScalars,
    const triangulationType *triangulation);

  template <typename scalarType, class triangulationType>
  int getPersistenceDiagram(
    vtkUnstructuredGrid *outputCTPersistenceDiagram,
    ttk::ftm::TreeType treeType,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram,
    vtkDataArray *inputScalars,
    const triangulationType *triangulation);

  template <typename scalarType, typename vtkSimplexArray, class triangulationType>
  int setPersistenceDiagramInfoInsideDomain(
    ttk::SimplexId id,
    vtkSmartPointer<vtkSimplexArray> vertexIdentifierScalars,
    vtkSmartPointer<vtkIntArray> nodeTypeScalars,
    vtkDataArray *birthScalars,
    vtkDataArray *deathScalars,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram,
    vtkSmartPointer<vtkPoints> points,
    vtkIdType ids[3],
    vtkDataArray *inputScalars,
    const triangulationType *triangulation);

  template <typename scalarType, class triangulationType>
  int getPersistenceDiagramInsideDomain(
    vtkUnstructuredGrid *outputCTPersistenceDiagram,
    ttk::ftm::TreeType treeType,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram,
    vtkDataArray *inputScalars,
    const triangulationType *triangulation);

  template <typename VTK_TT, typename TTK_TT>
  int dispatch(vtkUnstructuredGrid *outputCTPersistenceDiagram,
               vtkDataArray *inputScalarDataArray,
               const VTK_TT *inputScalars,
               int inputOffsetsDataType,
               const void *inputOffsets,
               const TTK_TT *triangulation);

  template <typename VTK_TT>
  int deleteDiagram();

protected:
  ttkPersistenceDiagram();
  ~ttkPersistenceDiagram() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  bool ForceInputOffsetScalarField{false};
  int ShowInsideDomain{false};

  bool computeDiagram_{true};
  void *CTDiagram_{nullptr};
  int scalarDataType{0};
};

#endif // _TTK_PERSISTENCEDIAGRAM_H
