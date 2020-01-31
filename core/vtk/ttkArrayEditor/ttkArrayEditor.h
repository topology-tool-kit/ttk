/// \ingroup vtk
/// \class ttkArrayEditor
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.10.2018
///
/// \brief TTK VTK-filter that adds/removes data to/from a vtkDataObject.
///
/// This filter adds/removes data arrays to/form a 'vtkDataObject' (called
/// target) based on a string and/or point/cell/field data of an optional second
/// 'vtkDataObject' (called source).

#pragma once

// VTK Module
#include <ttkArrayEditorModule.h>

// VTK includes
#include <ttkAlgorithm.h>

class TTKARRAYEDITOR_EXPORT ttkArrayEditor : public ttkAlgorithm {
private:
  int EditorMode{0};
  int TargetAttributeType{2};
  std::string DataString{""};
  std::vector<std::pair<int, std::string>> TargetArraySelection;

  std::string TargetArrayName;
  std::pair<int, std::string> TargetArray;
  int TargetArrayType;
  int TargetArrayIndexation[2];

  std::vector<std::pair<int, std::string>> SourceArraySelection;

public:
  static ttkArrayEditor *New();
  vtkTypeMacro(ttkArrayEditor, ttkAlgorithm);

  vtkSetMacro(EditorMode, int);
  vtkGetMacro(EditorMode, int);
  vtkSetMacro(TargetAttributeType, int);
  vtkGetMacro(TargetAttributeType, int);
  vtkSetMacro(DataString, std::string);
  vtkGetMacro(DataString, std::string);

  int SetTargetArray(
    int idx, int port, int connection, int arrayAssociation, const char *name) {
    this->TargetArray = {arrayAssociation, name};
    this->Modified();
    return 1;
  }
  vtkSetMacro(TargetArrayName, std::string);
  vtkGetMacro(TargetArrayName, std::string);
  vtkSetMacro(TargetArrayType, int);
  vtkGetMacro(TargetArrayType, int);
  vtkSetVector2Macro(TargetArrayIndexation, int);
  vtkGetVector2Macro(TargetArrayIndexation, int);

  /**
   * Adds an array to pass through.
   * fieldType where the array that should be passed (point data, cell data,
   * etc.). It should be one of the constants defined in the
   * vtkDataObject::AttributeTypes enumeration.
   */
  void
    AddArray(const int fieldType, const char *name, const int targetOrSource);
  void AddTargetPointDataArray(const char *name);
  void AddTargetCellDataArray(const char *name);
  void AddTargetFieldDataArray(const char *name);
  void AddSourcePointDataArray(const char *name);
  void AddSourceCellDataArray(const char *name);
  void AddSourceFieldDataArray(const char *name);

  void RemoveArray(const int fieldType,
                   const char *name,
                   bool deleteType,
                   const int targetOrSource);
  void RemoveTargetPointDataArray(const char *name);
  void RemoveTargetCellDataArray(const char *name);
  void RemoveTargetFieldDataArray(const char *name);
  void RemoveSourcePointDataArray(const char *name);
  void RemoveSourceCellDataArray(const char *name);
  void RemoveSourceFieldDataArray(const char *name);

  /**
   * Clear all arrays to pass through.
   */
  void ClearTargetPointDataArrays();
  void ClearTargetCellDataArrays();
  void ClearTargetFieldDataArrays();
  void ClearSourcePointDataArrays();
  void ClearSourceCellDataArrays();
  void ClearSourceFieldDataArrays();

protected:
  ttkArrayEditor();
  ~ttkArrayEditor();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};