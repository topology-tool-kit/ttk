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

  vtkSetMacro(InputScalarFieldName, std::string);
  vtkGetMacro(InputScalarFieldName, std::string);

  vtkSetMacro(InputIdentifiersFieldName, std::string);
  vtkGetMacro(InputIdentifiersFieldName, std::string);

  vtkSetMacro(OutputScalarFieldName, std::string);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(OutputScalarFieldType, int);
  vtkGetMacro(OutputScalarFieldType, int);

  vtkSetMacro(ForceInputScalarField, bool);
  vtkGetMacro(ForceInputScalarField, bool);

  vtkSetMacro(UseCotanWeights, bool);
  vtkGetMacro(UseCotanWeights, bool);

  vtkSetMacro(SolvingMethod, int);
  vtkGetMacro(SolvingMethod, int);

  // get mesh from VTK
  int getTriangulation(vtkDataSet *input);
  // get array of identifiers on the mesh
  int getIdentifiers(vtkPointSet *input);
  // get constraint values on identifiers
  int getConstraints(vtkPointSet *input);

  // default copy constructor
  ttkHarmonicFieldComputation(const ttkHarmonicFieldComputation &) = delete;
  // default move constructor
  ttkHarmonicFieldComputation(ttkHarmonicFieldComputation &&) = delete;
  // default copy assignment operator
  ttkHarmonicFieldComputation &operator=(const ttkHarmonicFieldComputation &)
    = delete;
  // default move assignment operator
  ttkHarmonicFieldComputation &operator=(ttkHarmonicFieldComputation &&)
    = delete;

protected:
  ttkHarmonicFieldComputation();

  ~ttkHarmonicFieldComputation() override = default;

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  // user-defined input constraints (float) scalar field name
  std::string InputScalarFieldName;
  // output scalar field
  std::string OutputScalarFieldName;
  // let the user choose a different identifier scalar field
  bool ForceInputScalarField;
  // graph laplacian variant
  bool UseCotanWeights;
  // user-defined input identifier (SimplexId) scalar field name
  std::string InputIdentifiersFieldName;
  // user-selected solving method
  int SolvingMethod;

  // enum: float or double
  int OutputScalarFieldType;
  // worker object
  ttk::HarmonicFieldComputation harmonicField_;
  // teh mesh
  ttk::Triangulation *triangulation_;
  // points on the mesh where constraints_ are set
  vtkDataArray *identifiers_;
  // scalar field constraint values on identifiers_
  vtkDataArray *constraints_;
};
