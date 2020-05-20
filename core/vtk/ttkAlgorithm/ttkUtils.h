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

namespace ttkUtils {
  int replaceVariable(const std::string &iString,
                      vtkFieldData *fieldData,
                      std::string &oString,
                      std::string &errorMsg);

  int replaceVariables(const std::string &iString,
                       vtkFieldData *fieldData,
                       std::string &oString,
                       std::string &errorMsg);

  int stringListToVector(const std::string &iString,
                         std::vector<std::string> &v);

  int stringListToDoubleVector(const std::string &iString,
                               std::vector<double> &v);

  vtkSmartPointer<vtkAbstractArray> csvToVtkArray(std::string line);

  vtkSmartPointer<vtkDoubleArray> csvToDoubleArray(std::string line);

  // Emultate old VTK functions

  void *GetVoidPointer(vtkDataArray *array, vtkIdType start = 0);
  void *GetVoidPointer(vtkPoints *points, vtkIdType start = 0);

  void *
    WriteVoidPointer(vtkDataArray *array, vtkIdType start, vtkIdType numValues);
  void *WritePointer(vtkDataArray *array, vtkIdType start, vtkIdType numValues);

  void SetVoidArray(vtkDataArray *array, void *data, vtkIdType size, int save);
}; // namespace ttkUtils