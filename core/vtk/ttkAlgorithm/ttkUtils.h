/// \ingroup vtk
/// \class ttkUtils
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.06.2019
///
/// \brief TTK Util Functions.

#pragma once

#include <vtkAbstractArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

class ttkUtils {
  public:

  ttkUtils();
  ~ttkUtils();

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

  static vtkSmartPointer<vtkAbstractArray> csvToVtkArray(std::string line);

  static vtkSmartPointer<vtkDoubleArray> csvToDoubleArray(std::string line);

  // Emultate old VTK functions

  static void *GetVoidPointer(vtkDataArray *array, vtkIdType start = 0);
  static void *GetVoidPointer(vtkPoints *points, vtkIdType start = 0);

  static void *
    WriteVoidPointer(vtkDataArray *array, vtkIdType start, vtkIdType numValues);
  static void *WritePointer(vtkDataArray *array, vtkIdType start, vtkIdType numValues);

  static void SetVoidArray(vtkDataArray *array, void *data, vtkIdType size, int save);
};