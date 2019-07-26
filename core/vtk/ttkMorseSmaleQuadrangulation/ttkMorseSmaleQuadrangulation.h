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

#pragma once

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>

// ttk code includes
#include <MorseSmaleQuadrangulation.h>
#include <ttkWrapper.h>

#include <ttkTriangulation.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkMorseSmaleQuadrangulation
#else
class ttkMorseSmaleQuadrangulation
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkMorseSmaleQuadrangulation *New();
  vtkTypeMacro(ttkMorseSmaleQuadrangulation, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkGetMacro(DualQuadrangulation, bool);
  vtkSetMacro(DualQuadrangulation, bool);

  vtkSetMacro(ShowResError, bool);
  vtkGetMacro(ShowResError, bool);

  // get data from Morse-Smale Complex
  int getCriticalPoints(vtkUnstructuredGrid *input);
  int getSeparatrices(vtkUnstructuredGrid *input);
  int getTriangulation(vtkUnstructuredGrid *input);

  // default copy constructor
  ttkMorseSmaleQuadrangulation(const ttkMorseSmaleQuadrangulation &) = delete;
  // default move constructor
  ttkMorseSmaleQuadrangulation(ttkMorseSmaleQuadrangulation &&) = delete;
  // default copy assignment operator
  ttkMorseSmaleQuadrangulation &operator=(const ttkMorseSmaleQuadrangulation &)
    = delete;
  // default move assignment operator
  ttkMorseSmaleQuadrangulation &operator=(ttkMorseSmaleQuadrangulation &&)
    = delete;

protected:
  ttkMorseSmaleQuadrangulation();

  ~ttkMorseSmaleQuadrangulation() override = default;

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  // triangulation
  ttk::Triangulation *triangulation_{};
  // if dual quadrangulation
  bool DualQuadrangulation{false};
  // show result despite error
  bool ShowResError{false};

  // worker object
  ttk::MorseSmaleQuadrangulation baseWorker_{};
};
