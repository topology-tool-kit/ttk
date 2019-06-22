#include <ttkJacobiSet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkJacobiSet)

  ttkJacobiSet::ttkJacobiSet() {
  ForceInputOffsetScalarField = false;

  // init
  UcomponentId = 0;
  VcomponentId = 1;

  UoffsetId = -1;
  VoffsetId = -1;

  EdgeIds = false;
  VertexScalars = false;
  UseAllCores = true;
}

ttkJacobiSet::~ttkJacobiSet() {
}

template <class dataTypeU, class dataTypeV>
int ttkJacobiSet::baseCall(vtkDataSet *input,
                           vtkDataArray *uField,
                           vtkDataArray *vField) {

  Timer t;

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;

  triangulation->setWrapper(this);
  JacobiSet<dataTypeU, dataTypeV> jacobiSet;
  jacobiSet.setWrapper(this);
  jacobiSet.setupTriangulation(triangulation);

  // point data
  jacobiSet.setInputField(uField->GetVoidPointer(0), vField->GetVoidPointer(0));

  vtkDataArray *offsetFieldU = NULL, *offsetFieldV = NULL;

  if((ForceInputOffsetScalarField)
     || ((UoffsetId != -1) && (VoffsetId != -1))) {
    if(OffsetFieldU.length()) {

      offsetFieldU = input->GetPointData()->GetArray(OffsetFieldU.data());

      if(offsetFieldU) {
        sosOffsetsU_.resize(offsetFieldU->GetNumberOfTuples());
        for(SimplexId i = 0; i < offsetFieldU->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU->GetTuple1(i);
        }

        jacobiSet.setSosOffsetsU(&sosOffsetsU_);
      }
    } else if(UoffsetId != -1) {
      offsetFieldU = input->GetPointData()->GetArray(UoffsetId);

      if(offsetFieldU) {
        sosOffsetsU_.resize(offsetFieldU->GetNumberOfTuples());
        for(SimplexId i = 0; i < offsetFieldU->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU->GetTuple1(i);
        }

        jacobiSet.setSosOffsetsU(&sosOffsetsU_);
      }
    } else if(input->GetPointData()->GetArray(ttk::OffsetFieldUName)) {
      offsetFieldU = input->GetPointData()->GetArray(ttk::OffsetFieldUName);

      if(offsetFieldU) {
        sosOffsetsU_.resize(offsetFieldU->GetNumberOfTuples());
        for(SimplexId i = 0; i < offsetFieldU->GetNumberOfTuples(); i++) {
          sosOffsetsU_[i] = offsetFieldU->GetTuple1(i);
        }

        jacobiSet.setSosOffsetsU(&sosOffsetsU_);
      }
    }
    if(OffsetFieldV.length()) {

      offsetFieldV = input->GetPointData()->GetArray(OffsetFieldV.data());

      if(offsetFieldV) {
        sosOffsetsV_.resize(offsetFieldV->GetNumberOfTuples());
        for(SimplexId i = 0; i < offsetFieldV->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV->GetTuple1(i);
        }

        jacobiSet.setSosOffsetsV(&sosOffsetsV_);
      }
    } else if(VoffsetId != -1) {
      offsetFieldV = input->GetPointData()->GetArray(VoffsetId);

      if(offsetFieldV) {
        sosOffsetsV_.resize(offsetFieldV->GetNumberOfTuples());
        for(SimplexId i = 0; i < offsetFieldV->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV->GetTuple1(i);
        }

        jacobiSet.setSosOffsetsV(&sosOffsetsV_);
      }
    } else if(input->GetPointData()->GetArray(ttk::OffsetFieldVName)) {
      offsetFieldV = input->GetPointData()->GetArray(ttk::OffsetFieldVName);

      if(offsetFieldV) {
        sosOffsetsV_.resize(offsetFieldV->GetNumberOfTuples());
        for(SimplexId i = 0; i < offsetFieldV->GetNumberOfTuples(); i++) {
          sosOffsetsV_[i] = offsetFieldV->GetTuple1(i);
        }

        jacobiSet.setSosOffsetsV(&sosOffsetsV_);
      }
    }
  }

  // go!
  jacobiSet.execute(jacobiSet_);
  Modified();

  return 0;
}
int ttkJacobiSet::doIt(vector<vtkDataSet *> &inputs,
                       vector<vtkDataSet *> &outputs) {

  Memory m;

  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  vtkDataArray *uComponent = NULL, *vComponent = NULL;

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

  {
    stringstream msg;
    msg << "[ttkJacobiSet] U-component: `" << uComponent->GetName() << "'"
        << endl;
    msg << "[ttkJacobiSet] V-component: `" << vComponent->GetName() << "'"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // set the jacobi functor
  switch(vtkTemplate2PackMacro(
    uComponent->GetDataType(), vComponent->GetDataType())) {
    ttkTemplate2Macro(
      baseCall<VTK_T1 TTK_COMMA VTK_T2>(input, uComponent, vComponent));
  }

  vtkSmartPointer<vtkCharArray> edgeTypes
    = vtkSmartPointer<vtkCharArray>::New();

  edgeTypes->SetNumberOfComponents(1);
  edgeTypes->SetNumberOfTuples(2 * jacobiSet_.size());
  edgeTypes->SetName("Critical Type");

  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  pointSet->SetNumberOfPoints(2 * jacobiSet_.size());

  vtkSmartPointer<vtkCellArray> cellArray
    = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(2);

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
  if(!triangulation)
    return -3;

  SimplexId pointCount = 0;
  double p[3];
  for(SimplexId i = 0; i < (SimplexId)jacobiSet_.size(); i++) {

    SimplexId edgeId = jacobiSet_[i].first;
    SimplexId vertexId0 = -1, vertexId1 = -1;
    triangulation->getEdgeVertex(edgeId, 0, vertexId0);
    triangulation->getEdgeVertex(edgeId, 1, vertexId1);

    input->GetPoint(vertexId0, p);
    pointSet->SetPoint(pointCount, p);
    edgeTypes->SetTuple1(pointCount, (float)jacobiSet_[i].second);
    idList->SetId(0, pointCount);
    pointCount++;

    input->GetPoint(vertexId1, p);
    pointSet->SetPoint(pointCount, p);
    edgeTypes->SetTuple1(pointCount, (float)jacobiSet_[i].second);
    idList->SetId(1, pointCount);
    pointCount++;

    cellArray->InsertNextCell(idList);
  }
  output->SetPoints(pointSet);
  output->SetCells(VTK_LINE, cellArray);
  output->GetPointData()->AddArray(edgeTypes);

  if(EdgeIds) {
    vtkSmartPointer<ttkSimplexIdTypeArray> edgeIdArray
      = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    edgeIdArray->SetNumberOfComponents(1);
    edgeIdArray->SetNumberOfTuples(jacobiSet_.size());
    edgeIdArray->SetName("EdgeIds");

    pointCount = 0;
    for(SimplexId i = 0; i < (SimplexId)jacobiSet_.size(); i++) {
      edgeIdArray->SetTuple1(pointCount, (float)jacobiSet_[i].first);
      pointCount++;
    }

    output->GetCellData()->AddArray(edgeIdArray);
  } else {
    output->GetCellData()->RemoveArray("EdgeIds");
  }

  if(VertexScalars) {

    for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {

      vtkDataArray *scalarField = input->GetPointData()->GetArray(i);

      switch(scalarField->GetDataType()) {

        case VTK_CHAR: {
          vtkSmartPointer<vtkCharArray> scalarArray
            = vtkSmartPointer<vtkCharArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(2 * jacobiSet_.size());
          scalarArray->SetName(scalarField->GetName());

          double *value = new double[scalarField->GetNumberOfComponents()];
          pointCount = 0;
          for(SimplexId j = 0; j < (SimplexId)jacobiSet_.size(); j++) {

            SimplexId edgeId = jacobiSet_[j].first;
            SimplexId vertexId0 = -1, vertexId1 = -1;
            triangulation->getEdgeVertex(edgeId, 0, vertexId0);
            triangulation->getEdgeVertex(edgeId, 1, vertexId1);

            scalarField->GetTuple(vertexId0, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;

            scalarField->GetTuple(vertexId1, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_DOUBLE: {
          vtkSmartPointer<vtkDoubleArray> scalarArray
            = vtkSmartPointer<vtkDoubleArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(2 * jacobiSet_.size());
          scalarArray->SetName(scalarField->GetName());

          double *value = new double[scalarField->GetNumberOfComponents()];
          pointCount = 0;
          for(SimplexId j = 0; j < (SimplexId)jacobiSet_.size(); j++) {

            SimplexId edgeId = jacobiSet_[j].first;
            SimplexId vertexId0 = -1, vertexId1 = -1;
            triangulation->getEdgeVertex(edgeId, 0, vertexId0);
            triangulation->getEdgeVertex(edgeId, 1, vertexId1);

            scalarField->GetTuple(vertexId0, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;

            scalarField->GetTuple(vertexId1, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_FLOAT: {
          vtkSmartPointer<vtkFloatArray> scalarArray
            = vtkSmartPointer<vtkFloatArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(2 * jacobiSet_.size());
          scalarArray->SetName(scalarField->GetName());

          double *value = new double[scalarField->GetNumberOfComponents()];
          pointCount = 0;
          for(SimplexId j = 0; j < (SimplexId)jacobiSet_.size(); j++) {

            SimplexId edgeId = jacobiSet_[j].first;
            SimplexId vertexId0 = -1, vertexId1 = -1;
            triangulation->getEdgeVertex(edgeId, 0, vertexId0);
            triangulation->getEdgeVertex(edgeId, 1, vertexId1);

            scalarField->GetTuple(vertexId0, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;

            scalarField->GetTuple(vertexId1, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_INT: {
          vtkSmartPointer<vtkIntArray> scalarArray
            = vtkSmartPointer<vtkIntArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(2 * jacobiSet_.size());
          scalarArray->SetName(scalarField->GetName());

          double *value = new double[scalarField->GetNumberOfComponents()];
          pointCount = 0;
          for(SimplexId j = 0; j < (SimplexId)jacobiSet_.size(); j++) {

            SimplexId edgeId = jacobiSet_[j].first;
            SimplexId vertexId0 = -1, vertexId1 = -1;
            triangulation->getEdgeVertex(edgeId, 0, vertexId0);
            triangulation->getEdgeVertex(edgeId, 1, vertexId1);

            scalarField->GetTuple(vertexId0, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;

            scalarField->GetTuple(vertexId1, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_ID_TYPE: {
          vtkSmartPointer<vtkIdTypeArray> scalarArray
            = vtkSmartPointer<vtkIdTypeArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(2 * jacobiSet_.size());
          scalarArray->SetName(scalarField->GetName());

          double *value = new double[scalarField->GetNumberOfComponents()];
          pointCount = 0;
          for(SimplexId j = 0; j < (SimplexId)jacobiSet_.size(); j++) {

            SimplexId edgeId = jacobiSet_[j].first;
            SimplexId vertexId0 = -1, vertexId1 = -1;
            triangulation->getEdgeVertex(edgeId, 0, vertexId0);
            triangulation->getEdgeVertex(edgeId, 1, vertexId1);

            scalarField->GetTuple(vertexId0, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;

            scalarField->GetTuple(vertexId1, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_UNSIGNED_SHORT: {
          vtkSmartPointer<vtkUnsignedShortArray> scalarArray
            = vtkSmartPointer<vtkUnsignedShortArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(2 * jacobiSet_.size());
          scalarArray->SetName(scalarField->GetName());

          double *value = new double[scalarField->GetNumberOfComponents()];
          pointCount = 0;
          for(SimplexId j = 0; j < (SimplexId)jacobiSet_.size(); j++) {

            SimplexId edgeId = jacobiSet_[j].first;
            SimplexId vertexId0 = -1, vertexId1 = -1;
            triangulation->getEdgeVertex(edgeId, 0, vertexId0);
            triangulation->getEdgeVertex(edgeId, 1, vertexId1);

            scalarField->GetTuple(vertexId0, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;

            scalarField->GetTuple(vertexId1, value);
            scalarArray->SetTuple(pointCount, value);
            pointCount++;
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        default: {
          stringstream msg;
          msg << "[ttkJacobiSet] Scalar attachment: "
              << "unsupported data type :(" << endl;
          dMsg(cerr, msg.str(), detailedInfoMsg);
        } break;
      }
    }
  } else {
    for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {
      output->GetPointData()->RemoveArray(
        input->GetPointData()->GetArray(i)->GetName());
    }
  }

  {
    stringstream msg;
    msg << "[ttkJacobiSet] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
