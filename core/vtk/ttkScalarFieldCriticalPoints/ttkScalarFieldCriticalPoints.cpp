#include <ttkScalarFieldCriticalPoints.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldCriticalPoints)

  ttkScalarFieldCriticalPoints::ttkScalarFieldCriticalPoints() {

  // init
  ForceInputOffsetScalarField = false;
  VertexBoundary = true;
  VertexIds = true;
  VertexScalars = true;

  ScalarFieldId = 0;
  OffsetFieldId = -1;
  OffsetField = ttk::OffsetScalarFieldName;

  UseAllCores = true;
}

ttkScalarFieldCriticalPoints::~ttkScalarFieldCriticalPoints() {
}

int ttkScalarFieldCriticalPoints::doIt(vector<vtkDataSet *> &inputs,
                                       vector<vtkDataSet *> &outputs) {
  Memory m;
  Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs.size()) {
    cerr
      << "[ttkScalarFieldCriticalPoints] Error: not enough input information."
      << endl;
    return -1;
  }
#endif

  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkScalarFieldCriticalPoints] Error: input pointer is NULL."
         << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkScalarFieldCriticalPoints] Error: input has no point." << endl;
    return -1;
  }
#endif

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    cerr << "[ttkScalarFieldCriticalPoints] Error: input triangulation is NULL."
         << endl;
    return -1;
  }
#endif

  if(VertexBoundary)
    triangulation->preprocessBoundaryVertices();

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField = NULL;
  vtkDataArray *offsetField = NULL;

  if(ScalarField.length()) {
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  } else {
    inputScalarField = input->GetPointData()->GetArray(ScalarFieldId);
  }

  if(!inputScalarField)
    return -1;

  {
    stringstream msg;
    msg << "[ttkScalarFieldCriticalPoints] Starting computation on field `"
        << inputScalarField->GetName() << "'..." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  if(OffsetFieldId != -1) {
    offsetField = input->GetPointData()->GetArray(OffsetFieldId);
    if(offsetField) {
      ForceInputOffsetScalarField = true;
      OffsetField = offsetField->GetName();
    }
  }

  if(ForceInputOffsetScalarField) {
    if(OffsetField.length()) {

      offsetField = input->GetPointData()->GetArray(OffsetField.data());
      // not good... in the future, we want to use the pointer itself...
      sosOffsets_.resize(offsetField->GetNumberOfTuples());
      for(SimplexId i = 0; i < offsetField->GetNumberOfTuples(); i++) {
        SimplexId offset = 0;
        offset = offsetField->GetTuple1(i);
        sosOffsets_[i] = offset;
      }
    }
  } else if(input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    offsetField = input->GetPointData()->GetArray(OffsetScalarFieldName);

    // not good... in the future, we want to use the pointer itself...
    sosOffsets_.resize(offsetField->GetNumberOfTuples());
    for(SimplexId i = 0; i < offsetField->GetNumberOfTuples(); i++) {
      SimplexId offset = 0;
      offset = offsetField->GetTuple1(i);
      sosOffsets_[i] = offset;
    }
  }

  switch(inputScalarField->GetDataType()) {

    vtkTemplateMacro({
      ScalarFieldCriticalPoints<VTK_TT> criticalPoints;
      criticalPoints.setupTriangulation(triangulation);
    });
  }

  int domainDimension = triangulation->getCellVertexNumber(0) - 1;

  switch(inputScalarField->GetDataType()) {
    vtkTemplateMacro({
      ScalarFieldCriticalPoints<VTK_TT> criticalPoints;

      criticalPoints.setWrapper(this);
      criticalPoints.setDomainDimension(domainDimension);
      // set up input
      // 1 -- vertex values
      criticalPoints.setScalarValues(inputScalarField->GetVoidPointer(0));
      criticalPoints.setVertexNumber(input->GetNumberOfPoints());

      // 2 -- set offsets (here, let the baseCode class fill it for us)
      criticalPoints.setSosOffsets(&sosOffsets_);

      // 3 -- set the connectivity
      criticalPoints.setupTriangulation(triangulation);

      // set up output
      criticalPoints.setOutput(&criticalPoints_);

      criticalPoints.execute();
    });
  }

  // allocate the output
  vtkSmartPointer<vtkCharArray> vertexTypes
    = vtkSmartPointer<vtkCharArray>::New();

  vertexTypes->SetNumberOfComponents(1);
  vertexTypes->SetNumberOfTuples(criticalPoints_.size());
  vertexTypes->SetName("CriticalType");

  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  pointSet->SetNumberOfPoints(criticalPoints_.size());
  double p[3];
  for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
    input->GetPoint(criticalPoints_[i].first, p);
    pointSet->SetPoint(i, p);
    vertexTypes->SetTuple1(i, (float)criticalPoints_[i].second);
  }
  output->SetPoints(pointSet);
  output->GetPointData()->AddArray(vertexTypes);

  if(VertexBoundary) {
    vtkSmartPointer<vtkCharArray> vertexBoundary
      = vtkSmartPointer<vtkCharArray>::New();
    vertexBoundary->SetNumberOfComponents(1);
    vertexBoundary->SetNumberOfTuples(criticalPoints_.size());
    vertexBoundary->SetName("IsOnBoundary");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
      vertexBoundary->SetTuple1(
        i, (char)triangulation->isVertexOnBoundary(criticalPoints_[i].first));
    }

    output->GetPointData()->AddArray(vertexBoundary);
  } else {
    output->GetPointData()->RemoveArray("IsOnBoundary");
  }

  if(VertexIds) {
    vtkSmartPointer<ttkSimplexIdTypeArray> vertexIds
      = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    vertexIds->SetNumberOfComponents(1);
    vertexIds->SetNumberOfTuples(criticalPoints_.size());
    vertexIds->SetName(ttk::VertexScalarFieldName);

    for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
      vertexIds->SetTuple1(i, criticalPoints_[i].first);
    }

    output->GetPointData()->AddArray(vertexIds);
  } else {
    output->GetPointData()->RemoveArray(ttk::VertexScalarFieldName);
  }

  if(VertexScalars) {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {

      vtkDataArray *scalarField = input->GetPointData()->GetArray(i);

      switch(scalarField->GetDataType()) {

        case VTK_CHAR: {
          vtkSmartPointer<vtkCharArray> scalarArray
            = vtkSmartPointer<vtkCharArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(criticalPoints_.size());
          scalarArray->SetName(scalarField->GetName());
          double *value = new double[scalarField->GetNumberOfComponents()];
          for(SimplexId j = 0; j < (SimplexId)criticalPoints_.size(); j++) {
            scalarField->GetTuple(criticalPoints_[j].first, value);
            scalarArray->SetTuple(j, value);
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_DOUBLE: {
          vtkSmartPointer<vtkDoubleArray> scalarArray
            = vtkSmartPointer<vtkDoubleArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(criticalPoints_.size());
          scalarArray->SetName(scalarField->GetName());
          double *value = new double[scalarField->GetNumberOfComponents()];
          for(SimplexId j = 0; j < (SimplexId)criticalPoints_.size(); j++) {
            scalarField->GetTuple(criticalPoints_[j].first, value);
            scalarArray->SetTuple(j, value);
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_FLOAT: {
          vtkSmartPointer<vtkFloatArray> scalarArray
            = vtkSmartPointer<vtkFloatArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(criticalPoints_.size());
          scalarArray->SetName(scalarField->GetName());
          double *value = new double[scalarField->GetNumberOfComponents()];
          for(SimplexId j = 0; j < (SimplexId)criticalPoints_.size(); j++) {
            scalarField->GetTuple(criticalPoints_[j].first, value);
            scalarArray->SetTuple(j, value);
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_INT: {
          vtkSmartPointer<vtkIntArray> scalarArray
            = vtkSmartPointer<vtkIntArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(criticalPoints_.size());
          scalarArray->SetName(scalarField->GetName());
          double *value = new double[scalarField->GetNumberOfComponents()];
          for(SimplexId j = 0; j < (SimplexId)criticalPoints_.size(); j++) {
            scalarField->GetTuple(criticalPoints_[j].first, value);
            scalarArray->SetTuple(j, value);
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        case VTK_ID_TYPE: {
          vtkSmartPointer<vtkIdTypeArray> scalarArray
            = vtkSmartPointer<vtkIdTypeArray>::New();
          scalarArray->SetNumberOfComponents(
            scalarField->GetNumberOfComponents());
          scalarArray->SetNumberOfTuples(criticalPoints_.size());
          scalarArray->SetName(scalarField->GetName());
          double *value = new double[scalarField->GetNumberOfComponents()];
          for(SimplexId j = 0; j < (SimplexId)criticalPoints_.size(); j++) {
            scalarField->GetTuple(criticalPoints_[j].first, value);
            scalarArray->SetTuple(j, value);
          }
          output->GetPointData()->AddArray(scalarArray);
          delete[] value;
        } break;

        default: {
          stringstream msg;
          msg << "[ttkScalarFieldCriticalPoints] Scalar attachment: "
              << "unsupported data type :(" << endl;
          dMsg(cerr, msg.str(), detailedInfoMsg);
        } break;
      }
    }
  } else {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {
      output->GetPointData()->RemoveArray(
        input->GetPointData()->GetArray(i)->GetName());
    }
  }

  {
    stringstream msg;
    msg << "[ttkScalarFieldCriticalPoints] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), 2);
  }

  return 0;
}
