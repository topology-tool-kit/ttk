/// \ingroup vtk
/// \class ttkArrayEditor
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.10.2018
///
/// \brief TTK VTK-filter that edit arrays of a vtkDataObject.
///
/// This filter adds data arrays to a 'vtkDataObject' (called
/// target) based on a string or point/cell/field data of an optional second
/// 'vtkDataObject' (called source). This filter can also be used to directly
/// edit an array (including renaming, type conversion, and reindexing).
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

#pragma once

// VTK Module
#include <ttkArrayEditorModule.h>

// VTK includes
#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>

class vtkDataArraySelection;

class TTKARRAYEDITOR_EXPORT ttkArrayEditor : public ttkAlgorithm {
public:
  enum class MODE {
    ADD_ARRAYS_FROM_STRING = 0,
    ADD_ARRAYS_FROM_SOURCE = 1,
    FILTER_ARRAYS_FROM_SOURCE = 2,
    EDIT_ARRAY = 3
  };

private:
  MODE EditorMode{MODE::ADD_ARRAYS_FROM_STRING};
  std::string DataString{""};
  bool ReplaceExistingArrays{true};

  std::string TargetArrayName;
  int TargetAssociation{2};
  int TargetArrayType;
  int TargetArrayIndexation[2];

  vtkSmartPointer<vtkDataArraySelection>
    ArraySelections[vtkDataObject::NUMBER_OF_ASSOCIATIONS];

public:
  static ttkArrayEditor *New();
  vtkTypeMacro(ttkArrayEditor, ttkAlgorithm);

  ttkSetEnumMacro(EditorMode, MODE);
  vtkGetEnumMacro(EditorMode, MODE);
  vtkSetMacro(TargetAssociation, int);
  vtkGetMacro(TargetAssociation, int);
  vtkSetMacro(DataString, const std::string &);
  vtkGetMacro(DataString, std::string);
  vtkSetMacro(ReplaceExistingArrays, bool);
  vtkGetMacro(ReplaceExistingArrays, bool);

  vtkSetMacro(TargetArrayName, const std::string &);
  vtkGetMacro(TargetArrayName, std::string);
  vtkSetMacro(TargetArrayType, int);
  vtkGetMacro(TargetArrayType, int);
  vtkSetVector2Macro(TargetArrayIndexation, int);
  vtkGetVector2Macro(TargetArrayIndexation, int);

  vtkDataArraySelection *GetArraySelection(int association);
  vtkDataArraySelection *GetPointDataArraySelection() {
    return this->GetArraySelection(vtkDataObject::FIELD_ASSOCIATION_POINTS);
  }
  vtkDataArraySelection *GetCellDataArraySelection() {
    return this->GetArraySelection(vtkDataObject::FIELD_ASSOCIATION_CELLS);
  }
  vtkDataArraySelection *GetFieldDataArraySelection() {
    return this->GetArraySelection(vtkDataObject::FIELD_ASSOCIATION_NONE);
  }
  vtkDataArraySelection *GetVertexDataArraySelection() {
    return this->GetArraySelection(vtkDataObject::FIELD_ASSOCIATION_VERTICES);
  }
  vtkDataArraySelection *GetEdgeDataArraySelection() {
    return this->GetArraySelection(vtkDataObject::FIELD_ASSOCIATION_EDGES);
  }
  vtkDataArraySelection *GetRowDataArraySelection() {
    return this->GetArraySelection(vtkDataObject::FIELD_ASSOCIATION_ROWS);
  }

protected:
  ttkArrayEditor();
  ~ttkArrayEditor();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
