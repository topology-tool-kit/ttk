#include <ttkArrayEditor.h>

#include <vtkInformation.h>

#include <vtkAbstractArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include <vtkCommand.h>
#include <vtkDataArraySelection.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkArrayEditor);

ttkArrayEditor::ttkArrayEditor() {
  this->setDebugMsgPrefix("ArrayEditor");

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  for(int cc = 0; cc < vtkDataObject::NUMBER_OF_ASSOCIATIONS; ++cc) {
    if(cc != vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS) {
      this->ArraySelections[cc] = vtkSmartPointer<vtkDataArraySelection>::New();
      this->ArraySelections[cc]->AddObserver(
        vtkCommand::ModifiedEvent, this, &ttkArrayEditor::Modified);
    } else {
      this->ArraySelections[cc] = nullptr;
    }
  }
}

ttkArrayEditor::~ttkArrayEditor() {
}

vtkDataArraySelection *ttkArrayEditor::GetArraySelection(int association) {
  if(association >= 0 && association < vtkDataObject::NUMBER_OF_ASSOCIATIONS) {
    return this->ArraySelections[association];
  }

  return nullptr;
}

int ttkArrayEditor::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    if(port == 1)
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkArrayEditor::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
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

  std::string associationNames[3] = {"point", "cell", "field"};

  // Pass Input to Output
  auto target = vtkDataObject::GetData(inputVector[0], 0);
  auto output = vtkDataObject::GetData(outputVector, 0);
  output->ShallowCopy(target);

  // switch based on mode
  if(this->EditorMode == 0) {
    this->printMsg(
      "Input string: '" + this->DataString + "'", ttk::debug::Priority::DETAIL);

    // From Data std::string
    if(this->DataString.length() > 0) {
      ttk::Timer t;

      // get target attribute
      int targetAssociation
        = this->TargetAssociation < 0 ? 2 : this->TargetAssociation;
      auto outputAtt = output->GetAttributesAsFieldData(targetAssociation);
      if(!outputAtt) {
        this->printErr("Target does not have requested attribute type.");
        return 0;
      }

      // if set to automatic then add to field data
      this->printMsg("Adding parsed arrays from string to "
                       + associationNames[targetAssociation] + " data",
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

      // Convert each line to a vtkDataArray
      {
        std::stringstream ss(finalExpressionString);
        std::string line;
        while(getline(ss, line, '\n')) {
          auto array = ttkUtils::csvToVtkArray(line);
          if(!array)
            return 0;

          if(this->ReplaceExistingArrays
             || !outputAtt->HasArray(array->GetName()))
            outputAtt->AddArray(array);
        }
      }

      this->printMsg("Adding parsed arrays from string to "
                       + associationNames[targetAssociation] + " data",
                     1, t.getElapsedTime());
    }
  } else if(this->EditorMode == 1) {
    ttk::Timer t;
    this->printMsg("Adding point/cell/field arrays from source", 0,
                   ttk::debug::LineMode::REPLACE);

    auto source = vtkDataObject::GetData(inputVector[1], 0);
    if(!source) {
      this->printErr("Unable to retrieve source vtkDataObject.");
      return 0;
    }

    for(int association = 0;
        association < vtkDataObject::NUMBER_OF_ASSOCIATIONS; ++association) {
      if(association == vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS)
        continue;

      // get source attributes and array selection
      auto sourceFD = source->GetAttributesAsFieldData(association);
      auto selection = this->GetArraySelection(association);
      if(!sourceFD || !selection)
        continue;

      // get target attribute
      int targetAssociation
        = this->TargetAssociation < 0 ? association : this->TargetAssociation;
      auto outputAtt = output->GetAttributesAsFieldData(targetAssociation);
      if(!outputAtt) {
        continue;
      }

      // add source arrays
      for(int i = 0; i < sourceFD->GetNumberOfArrays(); i++) {
        auto array = sourceFD->GetAbstractArray(i);

        if(!selection->ArrayIsEnabled(array->GetName()))
          continue;

        outputAtt->AddArray(array);
      }
    }

    this->printMsg(
      "Adding point/cell/field arrays from source", 1, t.getElapsedTime());

  } else if(this->EditorMode == 2) {
    ttk::Timer t;

    auto targetArray = this->GetInputArrayToProcess(0, inputVector);
    if(!targetArray) {
      this->printErr("Unable to retrieve input array.");
      return 0;
    }

    int targetArrayAssociation = this->GetInputArrayAssociation(0, inputVector);

    this->printMsg("Editing '" + std::string(targetArray->GetName()) + "' "
                     + associationNames[targetArrayAssociation] + " data array",
                   0, ttk::debug::LineMode::REPLACE);

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

    // get target attribute
    int targetAssociation = this->TargetAssociation < 0
                              ? targetArrayAssociation
                              : this->TargetAssociation;
    auto outputAtt = output->GetAttributesAsFieldData(targetAssociation);
    if(!outputAtt) {
      this->printErr("Target does not have requested attribute type.");
      return 0;
    }
    outputAtt->AddArray(copy);

    this->printMsg("Editing '" + std::string(targetArray->GetName()) + "' "
                     + associationNames[targetArrayAssociation] + " data array",
                   1, t.getElapsedTime());
  } else {
    this->printErr("Unsupported Editor Mode.");
    return 0;
  }

  return 1;
}