/// \ingroup vtk
/// \class ttkUtils
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.06.2019
///
/// \brief TTK Util Functions.

#pragma once

// VTK Module
#include <ttkAlgorithmModule.h>

#include <string>
#include <vector>

#include <vtkType.h>

class vtkFieldData;
class vtkDataArray;
class vtkDoubleArray;
class vtkPoints;
class vtkAbstractArray;
class vtkCellArray;

template <typename T>
class vtkSmartPointer;

class TTKALGORITHM_EXPORT ttkUtils {
private:
  ttkUtils() = default;

public:
  static int replaceVariable(const std::string &iString,
                             vtkFieldData *fieldData,
                             std::string &oString,
                             std::string &errorMsg);

  static int replaceVariables(const std::string &iString,
                              vtkFieldData *fieldData,
                              std::string &oString,
                              std::string &errorMsg);

  static int stringListToVector(const std::string &iString,
                                std::vector<std::string> &v);

  static int stringListToDoubleVector(const std::string &iString,
                                      std::vector<double> &v);

  static vtkSmartPointer<vtkAbstractArray>
    csvToVtkArray(const std::string &line);

  static vtkSmartPointer<vtkDoubleArray>
    csvToDoubleArray(const std::string &line);

  // Emultate old VTK functions
  static void *GetVoidPointer(vtkDataArray *array, vtkIdType start = 0);
  static void *GetVoidPointer(vtkPoints *points, vtkIdType start = 0);
  template <typename DT>
  static DT *GetPointer(vtkDataArray *array, vtkIdType start = 0) {
    return static_cast<DT *>(ttkUtils::GetVoidPointer(array, start));
  }

  static vtkSmartPointer<vtkAbstractArray> SliceArray(vtkAbstractArray *array,
                                                      vtkIdType idx);

  static void *
    WriteVoidPointer(vtkDataArray *array, vtkIdType start, vtkIdType numValues);
  static void *
    WritePointer(vtkDataArray *array, vtkIdType start, vtkIdType numValues);

  static void
    SetVoidArray(vtkDataArray *array, void *data, vtkIdType size, int save);

  // Fill Cell array using a pointer with the old memory layout
  // DEPRECTAED
  static void FillCellArrayFromSingle(vtkIdType const *cells,
                                      vtkIdType ncells,
                                      vtkCellArray *cellArray);

  // Fill Cell array using a pointer with the new memory layout
  static void FillCellArrayFromDual(vtkIdType const *cells_co,
                                    vtkIdType const *cells_off,
                                    vtkIdType ncells,
                                    vtkCellArray *cellArray);
};
