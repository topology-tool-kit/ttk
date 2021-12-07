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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///

#pragma once

// VTK Module
#include <ttkQuadrangulationSubdivisionModule.h>

// ttk code includes
#include <QuadrangulationSubdivision.h>
#include <ttkAlgorithm.h>

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

protected:
  ttkQuadrangulationSubdivision();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // display quadrangle statistics
  bool QuadStatistics{false};
};
