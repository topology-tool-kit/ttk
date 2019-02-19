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

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(OutputOffsetScalarFieldName, std::string);
  vtkGetMacro(OutputOffsetScalarFieldName, std::string);

  vtkSetMacro(InputVertexScalarFieldName, std::string);
  vtkGetMacro(InputVertexScalarFieldName, std::string);

  int getTriangulation(vtkDataSet *input);
  int getScalars(vtkDataSet *input);
  int getIdentifiers(vtkPointSet *input);
  int getOffsets(vtkDataSet *input);

protected:
  ttkHarmonicFieldComputation();

  ~ttkHarmonicFieldComputation();

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  int ScalarFieldId;
  int OffsetFieldId;
  std::string ScalarField;
  std::string InputOffsetScalarFieldName;
  std::string OutputOffsetScalarFieldName;
  std::string InputVertexScalarFieldName;
  bool hasUpdatedMesh_;

  ttk::HarmonicFieldComputation harmonicFieldComputation_;
  ttk::Triangulation *triangulation_;
  vtkDataArray *identifiers_;
  vtkDataArray *inputScalars_;
  vtkDataArray *offsets_;
  vtkDataArray *inputOffsets_;
};
