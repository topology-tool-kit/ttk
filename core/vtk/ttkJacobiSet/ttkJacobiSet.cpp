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
int ttkJacobiSet::dispatch(const dataTypeU *const uField,
                           const dataTypeV *const vField,
                           ttk::Triangulation *const triangulation) {
  ttkTemplateMacro(
    triangulation->getType(),
    this->execute(jacobiSet_, uField, vField,
                  *static_cast<TTK_TT *>(triangulation->getData()),
                  &isPareto_));
  return 0;
}

int ttkJacobiSet::RequestData(vtkInformation *ttkNotUsed(request),
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  const auto uComponent = this->GetInputArrayToProcess(0, input);
  const auto vComponent = this->GetInputArrayToProcess(1, input);

  if(uComponent == nullptr || vComponent == nullptr)
    return -1;

  this->printMsg("U-component: `" + std::string{uComponent->GetName()} + "'");
  this->printMsg("V-component: `" + std::string{vComponent->GetName()} + "'");

  // point data
  const auto offsetFieldU
    = this->GetOrderArray(input, 0, 2, ForceInputOffsetScalarField);
  const auto offsetFieldV
    = this->GetOrderArray(input, 1, 3, ForceInputOffsetScalarField);

  this->setSosOffsetsU(
    static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(offsetFieldU)));
  this->setSosOffsetsV(
    static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(offsetFieldV)));

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr)
    return -1;
  this->preconditionTriangulation(triangulation);

#ifndef TTK_ENABLE_DOUBLE_TEMPLATING
  if(uComponent->GetDataType() != vComponent->GetDataType()) {
    this->printErr(
      "Scalar fields should have same input type. Use TTKPointDataConverter or "
      "TTKArrayEditor to convert array types.");
    return 0;
  }
  switch(uComponent->GetDataType()) {
    vtkTemplateMacro(
      dispatch(static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(uComponent)),
               static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(vComponent)),
               triangulation));
  }
#else
  switch(vtkTemplate2PackMacro(
    uComponent->GetDataType(), vComponent->GetDataType())) {
    vtkTemplate2Macro(
      dispatch(static_cast<VTK_T1 *>(ttkUtils::GetVoidPointer(uComponent)),
               static_cast<VTK_T2 *>(ttkUtils::GetVoidPointer(vComponent)),
               triangulation));
  }
#endif // TTK_ENABLE_DOUBLE_TEMPLATING

  vtkNew<vtkSignedCharArray> edgeTypes{};

  edgeTypes->SetNumberOfComponents(1);
  edgeTypes->SetNumberOfTuples(2 * jacobiSet_.size());
  edgeTypes->SetName("Critical Type");

  vtkNew<vtkSignedCharArray> isPareto{};
  isPareto->SetNumberOfComponents(1);
  isPareto->SetNumberOfTuples(2 * jacobiSet_.size());
  isPareto->SetName("IsPareto");

  vtkNew<vtkPoints> pointSet{};
  pointSet->SetNumberOfPoints(2 * jacobiSet_.size());

  vtkNew<vtkCellArray> cellArray{};
  vtkNew<vtkIdList> idList{};
  idList->SetNumberOfIds(2);

  size_t pointCount = 0;
  std::array<double, 3> p{};
  for(size_t i = 0; i < jacobiSet_.size(); i++) {

    int edgeId = jacobiSet_[i].first;
    ttk::SimplexId vertexId0 = -1, vertexId1 = -1;
    triangulation->getEdgeVertex(edgeId, 0, vertexId0);
    triangulation->getEdgeVertex(edgeId, 1, vertexId1);

    input->GetPoint(vertexId0, p.data());
    pointSet->SetPoint(pointCount, p.data());
    edgeTypes->SetTuple1(pointCount, (float)jacobiSet_[i].second);
    isPareto->SetTuple1(pointCount, (float)isPareto_[i]);

    idList->SetId(0, pointCount);
    pointCount++;

    input->GetPoint(vertexId1, p.data());
    pointSet->SetPoint(pointCount, p.data());
    edgeTypes->SetTuple1(pointCount, (float)jacobiSet_[i].second);
    isPareto->SetTuple1(pointCount, (float)isPareto_[i]);
    idList->SetId(1, pointCount);
    pointCount++;

    cellArray->InsertNextCell(idList);
  }
  output->SetPoints(pointSet);
  output->SetCells(VTK_LINE, cellArray);
  output->GetPointData()->AddArray(edgeTypes);
  output->GetPointData()->AddArray(isPareto);

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
        ttk::SimplexId vertexId0 = -1, vertexId1 = -1;
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
