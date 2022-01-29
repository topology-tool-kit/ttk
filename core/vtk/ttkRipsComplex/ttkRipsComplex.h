/// \class ttkRipsComplex
/// \ingroup vtk
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date January 2022.
///
/// \brief TTK VTK-filter that wraps the ttk::RipsComplex processing
/// package.
///
/// VTK wrapping code for the ttk::RipsComplex package.
///
/// \param Input Input table (vtkTable)
/// \param Output Triangulation (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::RipsComplex

#pragma once

// VTK Module
#include <ttkRipsComplexModule.h>

// TTK includes
#include <RipsComplex.h>
#include <ttkAlgorithm.h>

class TTKRIPSCOMPLEX_EXPORT ttkRipsComplex : public ttkAlgorithm,
                                             protected ttk::RipsComplex {

public:
  static ttkRipsComplex *New();
  vtkTypeMacro(ttkRipsComplex, ttkAlgorithm);

  void SetScalarFields(const std::string &s) {
    ScalarFields.push_back(s);
    Modified();
  }

  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  vtkSetMacro(OutputDimension, int);
  vtkGetMacro(OutputDimension, int);

  vtkSetMacro(Epsilon, double);
  vtkGetMacro(Epsilon, double);

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(StdDev, double);
  vtkGetMacro(StdDev, double);

  vtkSetMacro(ComputeGaussianDensity, bool);
  vtkGetMacro(ComputeGaussianDensity, bool);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

  vtkSetMacro(XColumn, const std::string &);
  vtkGetMacro(XColumn, std::string);

  vtkSetMacro(YColumn, const std::string &);
  vtkGetMacro(YColumn, std::string);

  vtkSetMacro(ZColumn, const std::string &);
  vtkGetMacro(ZColumn, std::string);

protected:
  ttkRipsComplex();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool KeepAllDataArrays{true};
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> ScalarFields{};
  std::string XColumn{};
  std::string YColumn{};
  std::string ZColumn{};
};
