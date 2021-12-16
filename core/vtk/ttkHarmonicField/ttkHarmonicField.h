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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n

#pragma once

// VTK Module
#include <ttkHarmonicFieldModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
class vtkPointSet;

// ttk code includes
#include <HarmonicField.h>

class TTKHARMONICFIELD_EXPORT ttkHarmonicField : public ttkAlgorithm,
                                                 protected ttk::HarmonicField {

public:
  static ttkHarmonicField *New();

  vtkTypeMacro(ttkHarmonicField, ttkAlgorithm);

  vtkSetMacro(OutputScalarFieldName, const std::string &);
  vtkGetMacro(OutputScalarFieldName, std::string);

  vtkSetMacro(ForceConstraintIdentifiers, bool);
  vtkGetMacro(ForceConstraintIdentifiers, bool);

  vtkSetMacro(UseCotanWeights, bool);
  vtkGetMacro(UseCotanWeights, bool);

  void SetSolvingMethod(const int arg_) {
    if(arg_ == 0) {
      this->SolvingMethod = SolvingMethodUserType::AUTO;
    } else if(arg_ == 1) {
      this->SolvingMethod = SolvingMethodUserType::CHOLESKY;
    } else if(arg_ == 2) {
      this->SolvingMethod = SolvingMethodUserType::ITERATIVE;
    }
    this->Modified();
  }
  virtual int GetSolvingMethod() {
    switch(SolvingMethod) {
      case SolvingMethodUserType::AUTO:
        return 0;
      case SolvingMethodUserType::CHOLESKY:
        return 1;
      case SolvingMethodUserType::ITERATIVE:
        return 2;
    }
    return -1;
  }

  vtkSetMacro(LogAlpha, double);
  vtkGetMacro(LogAlpha, double);

protected:
  ttkHarmonicField();
  ~ttkHarmonicField() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // output scalar field
  std::string OutputScalarFieldName{"OutputHarmonicField"};
  // let the user choose a different identifier scalar field
  bool ForceConstraintIdentifiers{false};
  // graph laplacian variant
  bool UseCotanWeights{true};
  // user-selected solving method
  SolvingMethodUserType SolvingMethod{SolvingMethodUserType::AUTO};
  // penalty value
  double LogAlpha{5.0};

  // enum: float or double
  enum class FieldType { FLOAT, DOUBLE };
  FieldType OutputScalarFieldType{FieldType::FLOAT};
};
