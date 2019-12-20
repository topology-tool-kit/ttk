/// \ingroup vtk
/// \class ttkHarmonicField
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
/// \sa ttk::HarmonicField

#pragma once

// VTK includes -- to adapt
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

// VTK Module
#include <ttkHarmonicFieldModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// ttk code includes
#include <HarmonicField.h>

class TTKHARMONICFIELD_EXPORT ttkHarmonicField : public ttkAlgorithm,
                                                 protected ttk::HarmonicField {

public:
  static ttkHarmonicField *New();

  vtkTypeMacro(ttkHarmonicField, ttkAlgorithm);

  vtkSetMacro(InputScalarFieldName, std::string);
  vtkGetMacro(InputScalarFieldName, std::string);

  vtkSetMacro(InputIdentifiersFieldName, std::string);
  vtkGetMacro(InputIdentifiersFieldName, std::string);

  vtkSetMacro(OutputScalarFieldName, std::string);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(ForceConstraintIdentifiers, bool);
  vtkGetMacro(ForceConstraintIdentifiers, bool);

  vtkSetMacro(UseCotanWeights, bool);
  vtkGetMacro(UseCotanWeights, bool);

  vtkSetMacro(SolvingMethod, int);
  vtkGetMacro(SolvingMethod, int);

  vtkSetMacro(LogAlpha, double);
  vtkGetMacro(LogAlpha, double);

  // get array of identifiers on the mesh
  int getIdentifiers(vtkPointSet *input);
  // get constraint values on identifiers
  int getConstraints(vtkPointSet *input);

protected:
  ttkHarmonicField();

  ~ttkHarmonicField() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // user-defined input constraints (float) scalar field name
  std::string InputScalarFieldName{};
  // output scalar field
  std::string OutputScalarFieldName{"OutputHarmonicField"};
  // let the user choose a different identifier scalar field
  bool ForceConstraintIdentifiers{false};
  // graph laplacian variant
  bool UseCotanWeights{true};
  // user-defined input identifier (SimplexId) scalar field name
  std::string InputIdentifiersFieldName{ttk::VertexScalarFieldName};
  // user-selected solving method
  int SolvingMethod{0};
  // penalty value
  double LogAlpha{5};

  // enum: float or double
  enum class HarmonicFieldType { FLOAT, DOUBLE };
  HarmonicFieldType OutputScalarFieldType{HarmonicFieldType::FLOAT};

  // points on the mesh where constraints_ are set
  vtkDataArray *identifiers_{};
  // scalar field constraint values on identifiers_
  vtkDataArray *constraints_{};
};
