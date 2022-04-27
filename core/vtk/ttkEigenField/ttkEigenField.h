/// \ingroup vtk
/// \class ttkEigenField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date April 2019
///
/// \brief TTK VTK-filter for eigenfunctions computation.
///
/// This plugin computes the first eigenfunctions of a given
/// triangular surface mesh.
///
/// \param Input0 Input 2D geometry, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output Output eigenfunctions (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
///  "Spectral surface quadrangulation"
///  Shen Dong, Peer-Timo Bremer, Michael Garland, Valerio Pascucci, John C.
///  Hart SIGGRAPH 2006
///
/// \sa ttkHarmonicField
/// \sa ttk::EigenField
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///

#pragma once

// VTK Module
#include <ttkEigenFieldModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <EigenField.h>

class TTKEIGENFIELD_EXPORT ttkEigenField : public ttkAlgorithm,
                                           protected ttk::EigenField {

public:
  static ttkEigenField *New();

  vtkTypeMacro(ttkEigenField, ttkAlgorithm);

  vtkSetMacro(OutputFieldName, const std::string &);
  vtkGetMacro(OutputFieldName, std::string);

  vtkSetMacro(EigenNumber, unsigned int);
  vtkGetMacro(EigenNumber, unsigned int);

  vtkSetMacro(ComputeStatistics, bool);
  vtkGetMacro(ComputeStatistics, bool);

protected:
  ttkEigenField();
  ~ttkEigenField() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // output field name
  std::string OutputFieldName{"OutputEigenFunctions"};
  // number of eigenpairs to compute
  unsigned int EigenNumber{500};
  // if statistics are to be computed
  bool ComputeStatistics{false};

  // enum: float or double
  enum class FieldType { FLOAT, DOUBLE };

  FieldType OutputFieldType{FieldType::FLOAT};
};
