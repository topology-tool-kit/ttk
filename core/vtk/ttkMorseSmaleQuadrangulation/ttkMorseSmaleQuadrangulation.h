/// \ingroup vtk
/// \class ttkMorseSmaleQuadrangulation
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
/// \sa ttk::MorseSmaleQuadrangulation
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///

#pragma once

// VTK Module
#include <ttkMorseSmaleQuadrangulationModule.h>

// ttk code includes
#include <MorseSmaleQuadrangulation.h>
#include <ttkAlgorithm.h>

class TTKMORSESMALEQUADRANGULATION_EXPORT ttkMorseSmaleQuadrangulation
  : public ttkAlgorithm,
    virtual protected ttk::MorseSmaleQuadrangulation {

public:
  static ttkMorseSmaleQuadrangulation *New();
  vtkTypeMacro(ttkMorseSmaleQuadrangulation, ttkAlgorithm);

  vtkGetMacro(DualQuadrangulation, bool);
  vtkSetMacro(DualQuadrangulation, bool);

  vtkSetMacro(ShowResError, bool);
  vtkGetMacro(ShowResError, bool);

protected:
  ttkMorseSmaleQuadrangulation();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
