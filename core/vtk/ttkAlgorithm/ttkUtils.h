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

#include <limits>
#include <vtkStringArray.h>

int ttkUtils::replaceVariable(const std::string &iString,
                              vtkFieldData *fieldData,
                              std::string &oString,
                              std::string &errorMsg) {
  std::string varName = iString;
  int varIndex = -1;
  bool varIndexDefined = false;

  // Check if varIndex is specified
  size_t indexDelimiter0 = iString.find("[");
  size_t indexDelimiter1 = iString.find("]");
  if(indexDelimiter0 != std::string::npos
     && indexDelimiter1 != std::string::npos) {
    if(indexDelimiter0 > indexDelimiter1
       || iString.find("[", indexDelimiter0 + 1) != std::string::npos
       || iString.find("}", indexDelimiter1 + 1) != std::string::npos) {
      errorMsg = "Invalid Syntax:\n" + iString;
      return 0;
    }

    varName = iString.substr(0, indexDelimiter0);
    varIndex = stoi(iString.substr(
      indexDelimiter0 + 1, indexDelimiter1 - indexDelimiter0 - 1));
    varIndexDefined = true;
  }

  // Search in candidates
  auto column = fieldData->GetAbstractArray(varName.data());
  if(column == nullptr) {
    errorMsg = "FieldData does not contain array '" + varName + "'";
    return 0;
  }

  size_t n = column->GetNumberOfTuples();
  size_t m = column->GetNumberOfComponents();
  int s = n * m;

  if(!varIndexDefined) {
    if(s > 0) {
      oString = column->GetVariantValue(0).ToString();
      for(int i = 1; i < s; i++)
        oString += "," + column->GetVariantValue(i).ToString();
    }
  } else {
    if(varIndex < 0 || varIndex >= s) {
      errorMsg = "Index " + std::to_string(varIndex) + "/" + std::to_string(s)
                 + " for FieldData Array '" + varName + "' out of range";
      return 0;
    }
    oString = column->GetVariantValue(varIndex).ToString();
  }

  return 1;
}

int ttkUtils::replaceVariables(const std::string &iString,
                               vtkFieldData *fieldData,
                               std::string &oString,
                               std::string &errorMsg) {
  oString = iString;

  while(oString.find("{") != std::string::npos
        && oString.find("}") != std::string::npos) {
    size_t o = oString.find("{");
    size_t c = oString.find("}");
    // {...{....{...}...}..}
    // |            |
    // o            c

    size_t oNext = oString.find("{", o + 1);
    while(oNext != std::string::npos && oNext < c) {
      o = oNext;
      oNext = oString.find("{", o + 1);
    }

    // {...{....{var}...}..}
    //          |   |
    //          o   c
    std::string var = oString.substr(o + 1, c - 1 - o);

    std::string rVar;
    if(!replaceVariable(var, fieldData, rVar, errorMsg))
      return 0;

    oString = oString.substr(0, o) + rVar
              + oString.substr(c + 1, oString.length() - c - 1);
  }

  if(oString.find("{") != std::string::npos
     || oString.find("}") != std::string::npos) {
    errorMsg = "Invalid Syntax:\n" + iString;
    return 0;
  }

  return 1;
}

int ttkUtils::stringListToVector(const std::string &iString,
                                 std::vector<std::string> &v) {
  // ...,...,....,
  //     |  |
  //     i  j
  size_t i = 0;
  size_t j = iString.find(",");
  while(j != std::string::npos) {
    v.push_back(iString.substr(i, j - i));
    i = j + 1;
    j = iString.find(",", i);
  }
  if(iString.length() > i)
    v.push_back(iString.substr(i, iString.length() - i));

  return 1;
}

int ttkUtils::stringListToDoubleVector(const std::string &iString,
                                       std::vector<double> &v) {
  std::vector<std::string> stringVector;
  if(!ttkUtils::stringListToVector(iString, stringVector))
    return 0;

  size_t n = stringVector.size();
  v.resize(n);
  // try {
  for(size_t i = 0; i < n; i++)
    v[i] = stod(stringVector[i]);
  // } catch(std::invalid_argument &e) {
  //   return 0;
  // }

  return 1;
}

vtkSmartPointer<vtkAbstractArray> ttkUtils::csvToVtkArray(std::string line) {
  size_t firstComma = line.find(",", 0);

  if(firstComma == std::string::npos)
    return nullptr;

  std::string arrayName = line.substr(0, firstComma);

  std::vector<std::string> valuesAsString;
  stringListToVector(
    line.substr(firstComma + 1, std::string::npos), valuesAsString);
  size_t nValues = valuesAsString.size();
  if(nValues < 1)
    return nullptr;

  // Check if all elements are numbers
  bool isNumeric = true;
  {
    for(size_t i = 0; i < nValues; i++) {
      if(valuesAsString[i].size() > 0 && !isdigit(valuesAsString[i][0])) {
        isNumeric = false;
        break;
      }
    }
  }

  if(isNumeric) {
    auto array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetName(arrayName.data());
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(nValues);
    auto arrayData = reinterpret_cast<double *>(GetVoidPointer(array));
    // try {
    for(size_t i = 0; i < nValues; i++)
      arrayData[i] = std::stod(valuesAsString[i]);
    // } catch(std::invalid_argument &e) {
    //   // return nullptr;
    // }
    return array;
  } else {
    auto array = vtkSmartPointer<vtkStringArray>::New();
    array->SetName(arrayName.data());
    array->SetNumberOfValues(nValues);
    for(size_t i = 0; i < nValues; i++)
      array->SetValue(i, valuesAsString[i]);
    return array;
  }
}

vtkSmartPointer<vtkDoubleArray> ttkUtils::csvToDoubleArray(std::string line) {
  size_t firstComma = line.find(",", 0);

  if(firstComma == std::string::npos)
    return nullptr;

  std::string arrayName = line.substr(0, firstComma);
  std::string valuesAsString = line.substr(firstComma + 1, std::string::npos);

  std::vector<double> values;
  ttkUtils::stringListToDoubleVector(valuesAsString, values);
  size_t n = values.size();

  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(arrayName.data());
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(n);
  auto arrayData = reinterpret_cast<double *>(GetVoidPointer(array));
  for(size_t i = 0; i < n; i++)
    arrayData[i] = values[i];

  return array;
}

/// Retrieve pointer to the internal data
/// This method is a workaround to emulate
/// the old GetVoidPointer in vtkDataArray
void *ttkUtils::GetVoidPointer(vtkDataArray *array, vtkIdType start) {
  void *outPtr = nullptr;
  switch(array->GetDataType()) {
    vtkTemplateMacro(
      auto *aosArray = vtkAOSDataArrayTemplate<VTK_TT>::FastDownCast(array);
      if(aosArray) { outPtr = aosArray->GetVoidPointer(start); });
  }
  return outPtr;
}
void *ttkUtils::GetVoidPointer(vtkPoints *points, vtkIdType start) {
  return GetVoidPointer(points->GetData(), start);
}

void *ttkUtils::WriteVoidPointer(vtkDataArray *array,
                                 vtkIdType valueIdx,
                                 vtkIdType numValues) {
  void *outPtr = nullptr;
  switch(array->GetDataType()) {
    vtkTemplateMacro(auto *aosArray
                     = vtkAOSDataArrayTemplate<VTK_TT>::FastDownCast(array);
                     if(aosArray) {
                       outPtr = aosArray->WriteVoidPointer(valueIdx, numValues);
                     });
  }
  return outPtr;
}

void *ttkUtils::WritePointer(vtkDataArray *array,
                             vtkIdType valueIdx,
                             vtkIdType numValues) {
  void *outPtr = nullptr;
  switch(array->GetDataType()) {
    vtkTemplateMacro(
      auto *aosArray = vtkAOSDataArrayTemplate<VTK_TT>::FastDownCast(array);
      if(aosArray) { outPtr = aosArray->WritePointer(valueIdx, numValues); });
  }
  return outPtr;
}

void ttkUtils::SetVoidArray(vtkDataArray *array,
                            void *data,
                            vtkIdType size,
                            int save) {
  switch(array->GetDataType()) {
    vtkTemplateMacro(
      auto *aosArray = vtkAOSDataArrayTemplate<VTK_TT>::FastDownCast(array);
      if(aosArray) { aosArray->SetVoidArray(data, size, save); } else {
        std::cerr << "SetVoidArray on incompatible vtkDataArray:" << endl;
        array->Print(std::cerr);
      });
  }
}
