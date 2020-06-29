/// \ingroup vtk
/// \class ttkQuadrangulationSubdivision
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2019
///
/// \brief TTK VTK-filter for surface quadrangulation.
///
/// The current filter transforms a triangulated surface into a
/// quadrangulated one.
///
/// \param Input0 Input triangular surface (2D) geometry (vtkDataSet)
/// \param Output Quadrangular mesh (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkIntegralLines
/// \sa ttkFTMTree
/// \sa ttkIdentifiers
/// \sa ttk::QuadrangulationSubdivision

#pragma once

// VTK Module
#include <ttkQuadrangulationSubdivisionModule.h>

// ttk code includes
#include <QuadrangulationSubdivision.h>
#include <ttkAlgorithm.h>

class vtkUnstructuredGrid;

class TTKQUADRANGULATIONSUBDIVISION_EXPORT ttkQuadrangulationSubdivision
  : public ttkAlgorithm,
    protected ttk::QuadrangulationSubdivision {

public:
  static ttkQuadrangulationSubdivision *New();
  vtkTypeMacro(ttkQuadrangulationSubdivision, ttkAlgorithm);

  vtkSetMacro(SubdivisionLevel, unsigned int);
  vtkGetMacro(SubdivisionLevel, unsigned int);

  vtkSetMacro(RelaxationIterations, unsigned int);
  vtkGetMacro(RelaxationIterations, unsigned int);

  vtkSetMacro(LockInputExtrema, bool);
  vtkGetMacro(LockInputExtrema, bool);

  vtkSetMacro(LockAllInputVertices, bool);
  vtkGetMacro(LockAllInputVertices, bool);

  vtkSetMacro(ReverseProjection, bool);
  vtkGetMacro(ReverseProjection, bool);

  vtkSetMacro(HausdorffLevel, float);
  vtkGetMacro(HausdorffLevel, float);

  vtkSetMacro(ShowResError, bool);
  vtkGetMacro(ShowResError, bool);

  vtkSetMacro(QuadStatistics, bool);
  vtkGetMacro(QuadStatistics, bool);

  // default copy constructor
  ttkQuadrangulationSubdivision(const ttkQuadrangulationSubdivision &) = delete;
  // default move constructor
  ttkQuadrangulationSubdivision(ttkQuadrangulationSubdivision &&) = delete;
  // default copy assignment operator
  ttkQuadrangulationSubdivision &
    operator=(const ttkQuadrangulationSubdivision &)
    = delete;
  // default move assignment operator
  ttkQuadrangulationSubdivision &operator=(ttkQuadrangulationSubdivision &&)
    = delete;

protected:
  ttkQuadrangulationSubdivision();
  ~ttkQuadrangulationSubdivision() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int getQuadVertices(vtkUnstructuredGrid *input);

private:
  // number of subdivisions of the Morse-Smale Complex cells
  unsigned int SubdivisionLevel{1};
  // number of relaxation iterations
  unsigned int RelaxationIterations{10};
  // lock input extrema
  bool LockInputExtrema{false};
  // lock all input vertices
  bool LockAllInputVertices{false};
  // projection method
  bool ReverseProjection{false};
  // Hausdorff warning level
  float HausdorffLevel{200.F};
  // show result despite error
  bool ShowResError{false};
  // display quadrangle statistics
  bool QuadStatistics{false};
};
