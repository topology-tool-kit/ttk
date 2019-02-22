/// \ingroup vtk
/// \class ttkHarmonicFieldComputation
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date February 2019
///
/// \brief TTK VTK-filter for harmonic field computations.
///
/// The current filter takes a list of sources with attached scalar
/// values and produces a scalar harmonic field fulfilling these
/// constraints.
///
/// \param Input0 Input geometry, either 2D or 3D, either regular grid
/// or triangulation (vtkDataSet)
/// \param Input1 List of critical point constraints (vtkPointSet)
/// \param Output Output harmonic scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Dynamic harmonic fields for surface processing" \n
/// Kai Xu, Hao Zhang, Daniel Cohen-Or, Yueshan Xiong \n
/// Computers & Graphics 2009. \n
///
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkIntegralLines
/// \sa ttkFTMTree
/// \sa ttkIdentifiers
/// \sa ttk::HarmonicFieldComputation

#pragma once

// VTK includes -- to adapt
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
#include <HarmonicFieldComputation.h>
#include <ttkWrapper.h>

#include <ttkTriangulation.h>

enum HarmonicFieldType { Float = 0, Double };

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkHarmonicFieldComputation
#else
class ttkHarmonicFieldComputation
#endif
    : public vtkDataSetAlgorithm,
      public ttk::Wrapper {

public:
  static ttkHarmonicFieldComputation *New();
  vtkTypeMacro(ttkHarmonicFieldComputation, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(InputScalarFieldName, std::string);
  vtkGetMacro(InputScalarFieldName, std::string);

  vtkSetMacro(OutputScalarFieldName, std::string);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(OutputScalarFieldType, int);
  vtkGetMacro(OutputScalarFieldType, int);

  int getTriangulation(vtkDataSet *input);
  int getIdentifiers(vtkPointSet *input);

  // default copy constructor
  ttkHarmonicFieldComputation(const ttkHarmonicFieldComputation &) = delete;
  // default move constructor
  ttkHarmonicFieldComputation(ttkHarmonicFieldComputation &&) = delete;
  // default copy assignment operator
  ttkHarmonicFieldComputation &
  operator=(const ttkHarmonicFieldComputation &) = delete;
  // default move assignment operator
  ttkHarmonicFieldComputation &
  operator=(ttkHarmonicFieldComputation &&) = delete;

protected:
  ttkHarmonicFieldComputation();

  ~ttkHarmonicFieldComputation() override = default;

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  std::string ScalarField;
  std::string InputScalarFieldName;
  std::string OutputScalarFieldName;
  int OutputScalarFieldType;

  ttk::HarmonicFieldComputation harmonicField_;
  ttk::Triangulation *triangulation_;
  vtkDataArray *identifiers_;
};
