#include <ttkArrayEditor.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkInformationVector.h>

#include <vtkAbstractArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkArrayEditor);

ttkArrayEditor::ttkArrayEditor() {
  this->setDebugMsgPrefix("ArrayEditor");

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkArrayEditor::~ttkArrayEditor() {
}

int ttkArrayEditor::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    if(port == 1)
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  } else
    return 0;
  return 1;
}

int ttkArrayEditor::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

// =============================================================================
// GUI functions to add/remove/clear array selections
// =============================================================================

// Add Array Selections
void ttkArrayEditor::AddArray(const int fieldType,
                              const char *name,
                              const int targetOrSource) {
  std::string n(name);
  if(n.compare("") == 0) {
    this->printErr("Array name can not be ''");
    return;
  }

  if(targetOrSource == 0)
    this->TargetArraySelection.push_back(make_pair(fieldType, n));
  else
    this->SourceArraySelection.push_back(make_pair(fieldType, n));
  this->Modified();
}
void ttkArrayEditor::AddTargetPointDataArray(const char *name) {
  this->AddArray(vtkDataObject::POINT, name, 0);
}
void ttkArrayEditor::AddTargetCellDataArray(const char *name) {
  this->AddArray(vtkDataObject::CELL, name, 0);
}
void ttkArrayEditor::AddTargetFieldDataArray(const char *name) {
  this->AddArray(vtkDataObject::FIELD, name, 0);
}
void ttkArrayEditor::AddSourcePointDataArray(const char *name) {
  this->AddArray(vtkDataObject::POINT, name, 1);
}
void ttkArrayEditor::AddSourceCellDataArray(const char *name) {
  this->AddArray(vtkDataObject::CELL, name, 1);
}
void ttkArrayEditor::AddSourceFieldDataArray(const char *name) {
  this->AddArray(vtkDataObject::FIELD, name, 1);
}

// Remove Array Selections
void ttkArrayEditor::RemoveArray(const int fieldType,
                                 const char *name,
                                 bool deleteType,
                                 const int targetOrSource) {
  auto &selection = targetOrSource == 0 ? this->TargetArraySelection
                                        : this->SourceArraySelection;

  bool found = true;
  while(found) {
    found = false;
    auto iter = selection.begin();
    while(iter != selection.end()) {
      if(iter->first == fieldType && (deleteType || iter->second == name)) {
        found = true;
        iter = selection.erase(iter);
        this->Modified();
      } else
        ++iter;
    }
  }
}
void ttkArrayEditor::RemoveTargetPointDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::POINT, name, false, 0);
}
void ttkArrayEditor::RemoveTargetCellDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::CELL, name, false, 0);
}
void ttkArrayEditor::RemoveTargetFieldDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::FIELD, name, false, 0);
}
void ttkArrayEditor::RemoveSourcePointDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::POINT, name, false, 1);
}
void ttkArrayEditor::RemoveSourceCellDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::CELL, name, false, 1);
}
void ttkArrayEditor::RemoveSourceFieldDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::FIELD, name, false, 1);
}

// Clear Array Selections
void ttkArrayEditor::ClearTargetPointDataArrays() {
  this->RemoveArray(vtkDataObject::POINT, "", true, 0);
}
void ttkArrayEditor::ClearTargetCellDataArrays() {
  this->RemoveArray(vtkDataObject::CELL, "", true, 0);
}
void ttkArrayEditor::ClearTargetFieldDataArrays() {
  this->RemoveArray(vtkDataObject::FIELD, "", true, 0);
}
void ttkArrayEditor::ClearSourcePointDataArrays() {
  this->RemoveArray(vtkDataObject::POINT, "", true, 1);
}
void ttkArrayEditor::ClearSourceCellDataArrays() {
  this->RemoveArray(vtkDataObject::CELL, "", true, 1);
}
void ttkArrayEditor::ClearSourceFieldDataArrays() {
  this->RemoveArray(vtkDataObject::FIELD, "", true, 1);
}

template <typename VTK_T1, typename VTK_T2>
int copyArrayData(vtkDataArray *target, vtkDataArray *copy) {
  auto targetData = (VTK_T1 *)ttkUtils::GetVoidPointer(target);
  auto copyData = (VTK_T2 *)ttkUtils::GetVoidPointer(copy);
  for(size_t i = 0, n = target->GetNumberOfValues(); i < n; i++)
    copyData[i] = (VTK_T2)(targetData[i]);
  return 1;
}

// =============================================================================
// RequestData
// =============================================================================
int ttkArrayEditor::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {
  ttk::Timer globalTimer;

  std::string attributeTypeNames[3] = {"point", "cell", "field"};

  std::string modeString = this->EditorMode == 0
                             ? std::string("Add arrays parsed from std::string")
                             : this->EditorMode == 1
                                 ? std::string("Add arrays from source")
                                 : this->EditorMode == 2
                                     ? std::string("Remove arrays")
                                     : std::string("Edit Array");

  // Print Mode
  this->printMsg("Mode: " + modeString);
  this->printMsg(ttk::debug::Separator::L1);

  // Pass Input to Output
  auto target = vtkDataObject::GetData(inputVector[0], 0);
  auto source = vtkDataObject::GetData(inputVector[1], 0);

  auto output = vtkDataObject::GetData(outputVector, 0);

  output->ShallowCopy(target);
  auto outputAsDS = vtkDataSet::SafeDownCast(output);

  // Get relevant array selections
  auto &selection = this->EditorMode == 1 ? this->SourceArraySelection
                                          : this->TargetArraySelection;

  auto getarraysByAttributeType = [](vtkDataObject *object) {
    std::vector<vtkFieldData *> arraysByAttributeType(3, nullptr);

    auto objectAsDS = vtkDataSet::SafeDownCast(object);
    if(objectAsDS) {
      arraysByAttributeType[0] = (vtkFieldData *)objectAsDS->GetPointData();
      arraysByAttributeType[1] = (vtkFieldData *)objectAsDS->GetCellData();
    }
    arraysByAttributeType[2] = object->GetFieldData();

    return arraysByAttributeType;
  };

  auto checkData = [](std::vector<vtkFieldData *> &arraysByAttributeType,
                      std::string objectName, int attributeType,
                      vtkDataSet *objectAsDS, vtkAbstractArray *array,
                      ttkArrayEditor *caller) {
    if(!arraysByAttributeType[attributeType]) {
      caller->printErr(objectName + " does not have specified attribute type");
      return 0;
    }

    size_t n = attributeType == 0
                 ? objectAsDS->GetNumberOfPoints()
                 : attributeType == 1 ? objectAsDS->GetNumberOfCells() : 0;
    size_t m = array->GetNumberOfTuples();
    if(attributeType < 2 && n != m) {
      caller->printErr("tuple number of added array (" + std::to_string(m)
                       + ") != point/cell number of target ("
                       + std::to_string(n) + ")");
      return 0;
    }

    return 1;
  };

  // get output field type
  auto outputArraysByAttributeType = getarraysByAttributeType(output);

  // switch based on mode
  if(this->EditorMode == 0) {
    this->printMsg("Input std::string: '" + this->DataString + "'",
                   ttk::debug::Priority::DETAIL);

    // From Data std::string
    if(this->DataString.length() > 0) {
      ttk::Timer t;
      // if set to automatic then add to field data
      int targetAttributeType
        = this->TargetAttributeType < 0 ? 2 : this->TargetAttributeType;
      this->printMsg("Adding parsed arrays to "
                       + attributeTypeNames[targetAttributeType] + " data",
                     0, ttk::debug::LineMode::REPLACE);

      std::string finalExpressionString;
      {
        std::string errorMsg;
        if(!ttkUtils::replaceVariables(this->DataString, output->GetFieldData(),
                                       finalExpressionString, errorMsg)) {
          std::stringstream msg;
          this->printErr(errorMsg);
          return 0;
        }
      }

      {
        // For each line
        std::stringstream ss(finalExpressionString);
        std::string line;
        while(getline(ss, line, '\n')) {
          vtkSmartPointer<vtkAbstractArray> array
            = ttkUtils::csvToVtkArray(line);
          if(!array
             || !checkData(outputArraysByAttributeType, "Target",
                           targetAttributeType, outputAsDS, array, this))
            return 0;
          outputArraysByAttributeType[targetAttributeType]->AddArray(array);
        }
      }

      this->printMsg("Adding parsed arrays to "
                       + attributeTypeNames[targetAttributeType] + " data",
                     1, t.getElapsedTime());
    }
  } else if(this->EditorMode == 1) {
    ttk::Timer t;
    this->printMsg("Adding " + std::to_string(selection.size())
                     + " point/cell/field arrays from source",
                   0, ttk::debug::LineMode::REPLACE);

    if(inputVector[1]->GetNumberOfInformationObjects() == 1
       && selection.size() > 0) {

      auto source = inputVector[1]->GetInformationObject(0)->Get(
        vtkDataObject::DATA_OBJECT());
      auto sourcearraysByAttributeType = getarraysByAttributeType(source);

      for(size_t i = 0; i < selection.size(); i++) {
        auto &s = selection[i];

        // Get source array
        vtkFieldData *sourceFields = sourcearraysByAttributeType[s.first];
        vtkAbstractArray *sourceField
          = sourceFields ? vtkAbstractArray::SafeDownCast(
              sourceFields->GetAbstractArray(s.second.data()))
                         : nullptr;
        if(!sourceField) {
          this->printErr("Source does not have " + attributeTypeNames[s.first]
                         + " data array '" + s.second + "'");
          return 0;
        }

        // Copy source array
        auto copy
          = vtkSmartPointer<vtkAbstractArray>::Take(sourceField->NewInstance());
        auto copyAsDA = vtkDataArray::SafeDownCast(copy);
        auto copyAsSA = vtkStringArray::SafeDownCast(copy);

        if(copyAsDA)
          copyAsDA->ShallowCopy(vtkDataArray::SafeDownCast(sourceField));
        else if(copyAsSA) {
          copyAsSA->DeepCopy(vtkStringArray::SafeDownCast(sourceField));
          copyAsSA->SetName(sourceField->GetName());
        } else {
          this->printErr("Source array '" + s.second
                         + "' is not a 'vtkStringArray' nor a 'vtkDataArray'");
          return 0;
        }

        // Check if array can be safely added to target
        if(!checkData(outputArraysByAttributeType, "Target",
                      this->TargetAttributeType < 0 ? s.first
                                                    : this->TargetAttributeType,
                      outputAsDS, copy, this))
          return 0;
        outputArraysByAttributeType[this->TargetAttributeType < 0
                                      ? s.first
                                      : this->TargetAttributeType]
          ->AddArray(copy);
      }
    }

    this->printMsg("Adding " + std::to_string(selection.size())
                     + " point/cell/field arrays from source",
                   1, t.getElapsedTime());
  } else if(this->EditorMode == 2) {
    ttk::Timer t;
    this->printMsg("Removing " + std::to_string(selection.size())
                     + " point/cell/field arrays",
                   0, ttk::debug::LineMode::REPLACE);

    if(selection.size() > 0) {
      auto outputArraysByAttributeType = getarraysByAttributeType(output);

      for(size_t i = 0; i < selection.size(); i++) {
        auto &s = selection[i];

        vtkFieldData *outputFields = outputArraysByAttributeType[s.first];
        vtkAbstractArray *outputField
          = outputFields ? outputFields->GetAbstractArray(s.second.data())
                         : nullptr;
        if(!outputField) {
          this->printErr("Target does not have array '" + s.second
                         + "' of type '" + std::to_string(s.first) + "'");
          return 0;
        }

        outputFields->RemoveArray(s.second.data());
      }
    }

    this->printMsg("Removing " + std::to_string(selection.size())
                     + " point/cell/field arrays",
                   1, t.getElapsedTime());
  } else if(this->EditorMode == 3) {
    ttk::Timer t;
    this->printMsg("Editing '" + this->TargetArray.second + "' "
                     + attributeTypeNames[this->TargetArray.first]
                     + " data array",
                   0, ttk::debug::LineMode::REPLACE);

    auto targetArray = vtkDataArray::SafeDownCast(
      outputArraysByAttributeType[this->TargetArray.first]->GetAbstractArray(
        this->TargetArray.second.data()));
    if(!targetArray) {
      this->printErr("Only data arrays can be edited");
      return 0;
    }

    auto targetAttributeType = this->TargetAttributeType < 0
                                 ? this->TargetArray.first
                                 : this->TargetAttributeType;

    vtkSmartPointer<vtkDataArray> copy;
    // check if it is necessary to create a new array
    if(this->TargetArrayType >= 0 || this->TargetArrayIndexation[0] >= 0
       || this->TargetArrayIndexation[1] >= 0) {
      copy = vtkSmartPointer<vtkDataArray>::Take(vtkDataArray::SafeDownCast(
        vtkAbstractArray::CreateArray(this->TargetArrayType < 0
                                        ? targetArray->GetDataType()
                                        : this->TargetArrayType)));

      size_t nComponents = this->TargetArrayIndexation[1] < 0
                             ? targetArray->GetNumberOfComponents()
                             : this->TargetArrayIndexation[1];
      size_t nTuples = this->TargetArrayIndexation[0] < 0
                         ? targetArray->GetNumberOfTuples()
                         : this->TargetArrayIndexation[0];
      copy->Allocate(nTuples * nComponents);
      copy->SetNumberOfComponents(nComponents);
      copy->SetNumberOfTuples(nTuples);

      switch(vtkTemplate2PackMacro(
        targetArray->GetDataType(), copy->GetDataType())) {
        vtkTemplate2Macro((copyArrayData<VTK_T1, VTK_T2>(targetArray, copy)));
      }
    } else {
      copy = vtkSmartPointer<vtkDataArray>::Take(targetArray->NewInstance());
      copy->ShallowCopy(targetArray);
    }
    copy->SetName(this->TargetArrayName.compare("") == 0
                    ? targetArray->GetName()
                    : this->TargetArrayName.data());

    // Check if array can be safely added to target
    if(!checkData(outputArraysByAttributeType, "Target", targetAttributeType,
                  outputAsDS, copy, this))
      return 0;
    outputArraysByAttributeType[targetAttributeType]->AddArray(copy);

    this->printMsg("Editing '" + this->TargetArray.second + "' "
                     + attributeTypeNames[this->TargetArray.first]
                     + " data array",
                   1, t.getElapsedTime());
  }

  // Output Performance
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete", 1, globalTimer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
