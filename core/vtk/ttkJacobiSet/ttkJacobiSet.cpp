#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <ttkJacobiSet.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <array>

vtkStandardNewMacro(ttkJacobiSet);

ttkJacobiSet::ttkJacobiSet() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  ForceInputOffsetScalarField = false;
}

int ttkJacobiSet::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkJacobiSet::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <class dataTypeU, class dataTypeV>
int ttkJacobiSet::baseCall(vtkDataSet *input,
                           vtkDataArray *uField,
                           vtkDataArray *vField) {

  ttk::Timer t;

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(!triangulation)
    return -1;

  this->preconditionTriangulation(triangulation);

  // point data
  vtkDataArray *offsetFieldU = NULL, *offsetFieldV = NULL;

  if((ForceInputOffsetScalarField)
     || ((UoffsetId != -1) && (VoffsetId != -1))) {
    if(OffsetFieldU.length()) {

      offsetFieldU = input->GetPointData()->GetArray(OffsetFieldU.data());

      if(offsetFieldU) {
        sosOffsetsU_.resize(offsetFieldU->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldU->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU->GetTuple1(i);
        }

        this->setSosOffsetsU(&sosOffsetsU_);
      }
    } else if(UoffsetId != -1) {
      offsetFieldU = input->GetPointData()->GetArray(UoffsetId);

      if(offsetFieldU) {
        sosOffsetsU_.resize(offsetFieldU->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldU->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU->GetTuple1(i);
        }

        this->setSosOffsetsU(&sosOffsetsU_);
      }
    } else if(input->GetPointData()->GetArray(ttk::OffsetFieldUName)) {
      offsetFieldU = input->GetPointData()->GetArray(ttk::OffsetFieldUName);

      if(offsetFieldU) {
        sosOffsetsU_.resize(offsetFieldU->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldU->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU->GetTuple1(i);
        }

        this->setSosOffsetsU(&sosOffsetsU_);
      }
    }
    if(OffsetFieldV.length()) {

      offsetFieldV = input->GetPointData()->GetArray(OffsetFieldV.data());

      if(offsetFieldV) {
        sosOffsetsV_.resize(offsetFieldV->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldV->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV->GetTuple1(i);
        }

        this->setSosOffsetsV(&sosOffsetsV_);
      }
    } else if(VoffsetId != -1) {
      offsetFieldV = input->GetPointData()->GetArray(VoffsetId);

      if(offsetFieldV) {
        sosOffsetsV_.resize(offsetFieldV->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldV->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV->GetTuple1(i);
        }

        this->setSosOffsetsV(&sosOffsetsV_);
      }
    } else if(input->GetPointData()->GetArray(ttk::OffsetFieldVName)) {
      offsetFieldV = input->GetPointData()->GetArray(ttk::OffsetFieldVName);

      if(offsetFieldV) {
        sosOffsetsV_.resize(offsetFieldV->GetNumberOfTuples());
        for(vtkIdType i = 0; i < offsetFieldV->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV->GetTuple1(i);
        }

        this->setSosOffsetsV(&sosOffsetsV_);
      }
    }
  }

#define EXEC(_DATATYPE, TRIANGLCASE, TRIANGLTYPE, _CALL)                      \
  case TRIANGLCASE: {                                                         \
    this->execute(jacobiSet_,                                                 \
                  static_cast<dataTypeU *>(ttkUtils::GetVoidPointer(uField)), \
                  static_cast<dataTypeV *>(ttkUtils::GetVoidPointer(vField)), \
                  *static_cast<TRIANGLTYPE *>(triangulation->getData()));     \
  } break;

  ttkVtkTemplateTrianglMacro(_, triangulation->getType(), EXEC, _);

  Modified();

  return 0;
}

int ttkJacobiSet::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  vtkDataArray *uComponent = nullptr, *vComponent = nullptr;

  if(Ucomponent.length()) {
    uComponent = input->GetPointData()->GetArray(Ucomponent.data());
  } else {
    // default
    uComponent = input->GetPointData()->GetArray(UcomponentId);
  }
  if(!uComponent)
    return -1;

  if(Vcomponent.length()) {
    vComponent = input->GetPointData()->GetArray(Vcomponent.data());
  } else {
    // default
    vComponent = input->GetPointData()->GetArray(VcomponentId);
  }
  if(!vComponent)
    return -2;

  this->printMsg("U-component: `" + std::string{uComponent->GetName()} + "'");
  this->printMsg("V-component: `" + std::string{vComponent->GetName()} + "'");

  // set the jacobi functor
  switch(vtkTemplate2PackMacro(
    uComponent->GetDataType(), vComponent->GetDataType())) {
    vtkTemplate2Macro(
      (baseCall<VTK_T1, VTK_T2>(input, uComponent, vComponent)));
  }

  vtkNew<vtkSignedCharArray> edgeTypes{};

  edgeTypes->SetNumberOfComponents(1);
  edgeTypes->SetNumberOfTuples(2 * jacobiSet_.size());
  edgeTypes->SetName("Critical Type");

  vtkNew<vtkPoints> pointSet{};
  pointSet->SetNumberOfPoints(2 * jacobiSet_.size());

  vtkNew<vtkCellArray> cellArray{};
  vtkNew<vtkIdList> idList{};
  idList->SetNumberOfIds(2);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return -1;

  size_t pointCount = 0;
  std::array<double, 3> p{};
  for(size_t i = 0; i < jacobiSet_.size(); i++) {

    int edgeId = jacobiSet_[i].first;
    int vertexId0 = -1, vertexId1 = -1;
    triangulation->getEdgeVertex(edgeId, 0, vertexId0);
    triangulation->getEdgeVertex(edgeId, 1, vertexId1);

    input->GetPoint(vertexId0, p.data());
    pointSet->SetPoint(pointCount, p.data());
    edgeTypes->SetTuple1(pointCount, (float)jacobiSet_[i].second);
    idList->SetId(0, pointCount);
    pointCount++;

    input->GetPoint(vertexId1, p.data());
    pointSet->SetPoint(pointCount, p.data());
    edgeTypes->SetTuple1(pointCount, (float)jacobiSet_[i].second);
    idList->SetId(1, pointCount);
    pointCount++;

    cellArray->InsertNextCell(idList);
  }
  output->SetPoints(pointSet);
  output->SetCells(VTK_LINE, cellArray);
  output->GetPointData()->AddArray(edgeTypes);

  if(EdgeIds) {
    vtkNew<ttkSimplexIdTypeArray> edgeIdArray{};
    edgeIdArray->SetNumberOfComponents(1);
    edgeIdArray->SetNumberOfTuples(jacobiSet_.size());
    edgeIdArray->SetName("EdgeIds");

    pointCount = 0;
    for(size_t i = 0; i < jacobiSet_.size(); i++) {
      edgeIdArray->SetTuple1(pointCount, (float)jacobiSet_[i].first);
      pointCount++;
    }

    output->GetCellData()->AddArray(edgeIdArray);
  } else {
    output->GetCellData()->RemoveArray("EdgeIds");
  }

  if(VertexScalars) {

    for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {

      const auto scalarField = input->GetPointData()->GetArray(i);
      vtkSmartPointer<vtkDataArray> scalarArray{scalarField->NewInstance()};

      scalarArray->SetNumberOfComponents(scalarField->GetNumberOfComponents());
      scalarArray->SetNumberOfTuples(2 * jacobiSet_.size());
      scalarArray->SetName(scalarField->GetName());
      std::vector<double> value(scalarField->GetNumberOfComponents());

      for(size_t j = 0; j < jacobiSet_.size(); j++) {
        int edgeId = jacobiSet_[j].first;
        int vertexId0 = -1, vertexId1 = -1;
        triangulation->getEdgeVertex(edgeId, 0, vertexId0);
        triangulation->getEdgeVertex(edgeId, 1, vertexId1);

        scalarField->GetTuple(vertexId0, value.data());
        scalarArray->SetTuple(2 * j, value.data());

        scalarField->GetTuple(vertexId1, value.data());
        scalarArray->SetTuple(2 * j + 1, value.data());
      }
      output->GetPointData()->AddArray(scalarArray);
    }
  } else {
    for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {
      output->GetPointData()->RemoveArray(
        input->GetPointData()->GetArray(i)->GetName());
    }
  }

  return 1;
}
