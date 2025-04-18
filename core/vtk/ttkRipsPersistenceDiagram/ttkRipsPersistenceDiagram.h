/// \ingroup vtk
/// \class ttkRipsPersistenceDiagram
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.
///
/// \brief TTK VTK-filter that wraps the ttk::RipsPersistenceDiagram module.
///
/// VTK wrapping code for the ttk::RipsPersistenceDiagram package.
///
/// \param Input Input table (vtkTable)
/// \param Output PersistenceDiagram (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::RipsPersistenceDiagram
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkRipsPersistenceDiagramModule.h>

// VTK Includes
#include <ttkMacros.h>
#include <vtkUnstructuredGrid.h>

// TTK Includes
#include <RipsPersistenceDiagram.h>
#include <ttkAlgorithm.h>

/**
 * @brief Converts a Rips Persistence Diagram in the
 * ttk::rpd::MultidimensionalDiagram format to the VTK Unstructured Grid format.
 *
 * @param[out] vtu Output VTK Unstructured Grid
 * @param[in] diagram ttk::rpd::MultidimensionalDiagram to be converted
 * @param[in] SimplexMaximumDiameter Maximum diameter of simplices (to cap
 * infinite pairs)
 */
TTKRIPSPERSISTENCEDIAGRAM_EXPORT void
  DiagramToVTU(vtkUnstructuredGrid *vtu,
               const ttk::rpd::MultidimensionalDiagram &diagram,
               double SimplexMaximumDiameter);

/**
 * @brief Converts a vector of 1-dimensional persistent generators in the
 * ttk::rpd::Generator format to the VTK Unstructured Grid format.
 *
 * @param[out] vtu Output VTK Unstructured Grid
 * @param[in] inputPoints vtkPoints used to embed the generators
 * @param[in] generators Vector of ttk::rpd::Generator to be converted
 * @param[in] parametrize Whether to parametrize the generators
 */
TTKRIPSPERSISTENCEDIAGRAM_EXPORT void
  GeneratorsToVTU(vtkUnstructuredGrid *vtu,
                  vtkPoints *inputPoints,
                  const std::vector<ttk::rpd::Generator> &generators,
                  bool parametrize = true);

class TTKRIPSPERSISTENCEDIAGRAM_EXPORT ttkRipsPersistenceDiagram
  : public ttkAlgorithm, // we inherit from the generic ttkAlgorithm class
    protected ttk::RipsPersistenceDiagram { // and we inherit from the base
                                            // class
private:
  bool KeepAllDataArrays{true};
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> ScalarFields{};

public:
  static ttkRipsPersistenceDiagram *New();
  vtkTypeMacro(ttkRipsPersistenceDiagram, ttkAlgorithm);

  void SetScalarFields(const std::string &s) {
    ScalarFields.push_back(s);
    Modified();
  }

  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

  ttkSetEnumMacro(BackEnd, BACKEND);
  vtkGetEnumMacro(BackEnd, BACKEND);

  vtkSetMacro(SimplexMaximumDimension, int);
  vtkGetMacro(SimplexMaximumDimension, int);

  vtkSetMacro(SimplexMaximumDiameter, double);
  vtkGetMacro(SimplexMaximumDiameter, double);

  vtkSetMacro(FieldOfCoefficients, int);
  vtkGetMacro(FieldOfCoefficients, int);

  vtkSetMacro(InputIsDistanceMatrix, bool);
  vtkGetMacro(InputIsDistanceMatrix, bool);

  vtkSetMacro(DelaunayRips, bool);
  vtkGetMacro(DelaunayRips, bool);

  vtkSetMacro(OutputGenerators, bool);
  vtkGetMacro(OutputGenerators, bool);

protected:
  ttkRipsPersistenceDiagram();
  ~ttkRipsPersistenceDiagram() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};