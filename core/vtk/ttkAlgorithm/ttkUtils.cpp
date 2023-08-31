#include <ttkUtils.h>

#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

#ifdef TTK_ENABLE_MPI_TIME
#include <mpi.h>
#endif
int ttkUtils::replaceVariable(const std::string &iString,
                              vtkFieldData *fieldData,
                              std::string &oString,
                              std::string &errorMsg) {
  std::string varName = iString;
  int varIndex = -1;
  bool varIndexDefined = false;

  // Check if varIndex is specified
  size_t const indexDelimiter0 = iString.find('[');
  size_t const indexDelimiter1 = iString.find(']');
  if(indexDelimiter0 != std::string::npos
     && indexDelimiter1 != std::string::npos) {
    if(indexDelimiter0 > indexDelimiter1
       || iString.find('[', indexDelimiter0 + 1) != std::string::npos
       || iString.find('}', indexDelimiter1 + 1) != std::string::npos) {
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

  size_t const n = column->GetNumberOfTuples();
  size_t const m = column->GetNumberOfComponents();
  int const s = n * m;

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

  while(oString.find('{') != std::string::npos
        && oString.find('}') != std::string::npos) {
    size_t o = oString.find('{');
    size_t const c = oString.find('}');
    // {...{....{...}...}..}
    // |            |
    // o            c

    size_t oNext = oString.find('{', o + 1);
    while(oNext != std::string::npos && oNext < c) {
      o = oNext;
      oNext = oString.find('{', o + 1);
    }

    // {...{....{var}...}..}
    //          |   |
    //          o   c
    std::string const var = oString.substr(o + 1, c - 1 - o);

    std::string rVar;
    if(!replaceVariable(var, fieldData, rVar, errorMsg))
      return 0;

    oString = oString.substr(0, o).append(rVar).append(
      oString.substr(c + 1, oString.length() - c - 1));
  }

  if(oString.find('{') != std::string::npos
     || oString.find('}') != std::string::npos) {
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
  size_t j = iString.find(',');
  while(j != std::string::npos) {
    v.push_back(iString.substr(i, j - i));
    i = j + 1;
    j = iString.find(',', i);
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

  size_t const n = stringVector.size();
  v.resize(n);
  // try {
  for(size_t i = 0; i < n; i++)
    v[i] = stod(stringVector[i]);
  // } catch(std::invalid_argument &e) {
  //   return 0;
  // }

  return 1;
}

vtkSmartPointer<vtkAbstractArray>
  ttkUtils::csvToVtkArray(const std::string &line) {
  size_t const firstComma = line.find(',', 0);

  if(firstComma == std::string::npos)
    return nullptr;

  std::string const arrayName = line.substr(0, firstComma);

  std::vector<std::string> valuesAsString;
  stringListToVector(
    line.substr(firstComma + 1, std::string::npos), valuesAsString);
  size_t const nValues = valuesAsString.size();
  if(nValues < 1)
    return nullptr;

  // Check if all elements are numbers
  bool isNumeric = true;

  std::vector<double> valuesAsDouble(nValues);
  try {
    for(size_t i = 0; i < nValues; i++)
      valuesAsDouble[i] = std::stod(valuesAsString[i]);
  } catch(const std::invalid_argument &) {
    isNumeric = false;
  } catch(const std::out_of_range &) {
    isNumeric = false;
  }

  if(isNumeric) {
    auto array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetName(arrayName.data());
    array->SetNumberOfValues(nValues);
    for(size_t i = 0; i < nValues; i++)
      array->SetValue(i, valuesAsDouble[i]);
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

vtkSmartPointer<vtkDoubleArray>
  ttkUtils::csvToDoubleArray(const std::string &line) {
  size_t const firstComma = line.find(',', 0);

  if(firstComma == std::string::npos)
    return nullptr;

  std::string const arrayName = line.substr(0, firstComma);
  std::string const valuesAsString
    = line.substr(firstComma + 1, std::string::npos);

  std::vector<double> values;
  ttkUtils::stringListToDoubleVector(valuesAsString, values);
  size_t const n = values.size();

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
  if(array == nullptr)
    return outPtr;

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

vtkSmartPointer<vtkAbstractArray> ttkUtils::SliceArray(vtkAbstractArray *array,
                                                       vtkIdType idx) {
  auto slicedArray
    = vtkSmartPointer<vtkAbstractArray>::Take(array->NewInstance());
  slicedArray->SetName(array->GetName());
  slicedArray->SetNumberOfComponents(array->GetNumberOfComponents());
  slicedArray->SetNumberOfTuples(1);
  slicedArray->SetTuple(0, idx, array);
  return slicedArray;
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

[[deprecated]] void ttkUtils::FillCellArrayFromSingle(vtkIdType const *cells,
                                                      vtkIdType ncells,
                                                      vtkCellArray *cellArray) {
  size_t curPos = 0;
  vtkNew<vtkIdList> verts;
  for(vtkIdType cid = 0; cid < ncells; cid++) {
    const vtkIdType nbVerts = cells[curPos];
    verts->SetNumberOfIds(nbVerts);
    curPos++;
    for(vtkIdType v = 0; v < nbVerts; v++) {
      verts->SetId(v, cells[curPos]);
      curPos++;
    }
    cellArray->InsertNextCell(verts);
  }
}

void ttkUtils::FillCellArrayFromDual(vtkIdType const *cells_co,
                                     vtkIdType const *cells_off,
                                     vtkIdType ncells,
                                     vtkCellArray *cellArray) {
  size_t curPos = 0;
  vtkNew<vtkIdList> verts;
  for(vtkIdType cid = 0; cid < ncells; cid++) {
    const vtkIdType nbVerts = cells_off[cid + 1] - cells_off[cid];
    verts->SetNumberOfIds(nbVerts);
    for(vtkIdType v = 0; v < nbVerts; v++) {
      verts->SetId(v, cells_co[curPos]);
      curPos++;
    }
    cellArray->InsertNextCell(verts);
  }
}

int ttkUtils::CellVertexFromPoints(vtkDataSet *const dataSet,
                                   vtkPoints *const points) {

  if(dataSet == nullptr || points == nullptr) {
    return 0;
  }

  if(!dataSet->IsA("vtkUnstructuredGrid") && !dataSet->IsA("vtkPolyData")) {
    return 0;
  }

  const size_t nPoints = points->GetNumberOfPoints();
  if(nPoints == 0) {
    return 0;
  }

  vtkNew<vtkIdTypeArray> offsets{};
  offsets->SetNumberOfTuples(nPoints + 1);
  auto offsetsData = ttkUtils::GetPointer<vtkIdType>(offsets);
  for(size_t i = 0; i <= nPoints; i++)
    offsetsData[i] = i;

  vtkNew<vtkIdTypeArray> connectivity{};
  connectivity->SetNumberOfTuples(nPoints);
  auto connectivityData = ttkUtils::GetPointer<vtkIdType>(connectivity);
  for(size_t i = 0; i < nPoints; i++)
    connectivityData[i] = i;

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);

  if(dataSet->IsA("vtkUnstructuredGrid")) {
    const auto vtu{vtkUnstructuredGrid::SafeDownCast(dataSet)};
    if(vtu != nullptr) {
      vtu->SetPoints(points);
      vtu->SetCells(VTK_VERTEX, cells);
    }
  } else if(dataSet->IsA("vtkPolyData")) {
    const auto vtp{vtkPolyData::SafeDownCast(dataSet)};
    if(vtp != nullptr) {
      vtp->SetPoints(points);
      vtp->SetVerts(cells);
    }
  }

  return 1;
}
