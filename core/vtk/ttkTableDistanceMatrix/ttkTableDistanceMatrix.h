/// \ingroup vtk
/// \class ttkTableDistanceMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date January 2022
///
/// \brief Computes a distance matrix using LDistance from a vtkTable
///
/// \sa LDistanceMatrix
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_householdAnalysis/">Persistent
///   Generators Household Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_periodicPicture/">Persistent
///   Generators Periodic Picture example</a> \n

#pragma once

// VTK Module
#include <ttkTableDistanceMatrixModule.h>

// TTK code includes
#include <LDistanceMatrix.h>
#include <ttkAlgorithm.h>

class TTKTABLEDISTANCEMATRIX_EXPORT ttkTableDistanceMatrix
  : public ttkAlgorithm,
    protected ttk::LDistanceMatrix {

public:
  static ttkTableDistanceMatrix *New();

  vtkTypeMacro(ttkTableDistanceMatrix, ttkAlgorithm);

  void SetScalarFields(const std::string &s) {
    ScalarFields.emplace_back(s);
    Modified();
  }

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  vtkSetMacro(DistanceType, const std::string &);
  vtkGetMacro(DistanceType, std::string);

protected:
  ttkTableDistanceMatrix();
  ~ttkTableDistanceMatrix() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::vector<std::string> ScalarFields{};
  std::string RegexpString{".*"};
  bool SelectFieldsWithRegexp{false};
};
