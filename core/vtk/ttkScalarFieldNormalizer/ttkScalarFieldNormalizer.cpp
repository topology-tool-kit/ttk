#include <ttkScalarFieldNormalizer.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldNormalizer)

  ttkScalarFieldNormalizer::ttkScalarFieldNormalizer() {

  // init
  outputScalarField_ = NULL;

  UseAllCores = true;
}

ttkScalarFieldNormalizer::~ttkScalarFieldNormalizer() {

  if(outputScalarField_)
    outputScalarField_->Delete();
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkScalarFieldNormalizer::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkScalarFieldNormalizer::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkScalarFieldNormalizer] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkScalarFieldNormalizer::normalize(vtkDataArray *input,
                                        vtkDataArray *output) const {

  if(!output)
    return -1;
  if(!input)
    return -2;

  double min = 0, max = 0;
  for(SimplexId i = 0; i < input->GetNumberOfTuples(); i++) {

    double value = input->GetTuple1(i);

    if((!i) || (value < min)) {
      min = value;
    }
    if((!i) || (value > max)) {
      max = value;
    }
  }

  for(SimplexId i = 0; i < input->GetNumberOfTuples(); i++) {
    double value = input->GetTuple1(i);

    value = (value - min) / (max - min) + pow10(-FLT_DIG);

    output->SetTuple1(i, value);
  }

  return 0;
}

int ttkScalarFieldNormalizer::doIt(vtkDataSet *input, vtkDataSet *output) {

  // test
  //   ttkTriangulation myTriangulation;
  //   myTriangulation.setInputData(input);
  //
  //   myTriangulation.preprocessCellTriangles();
  //
  //   for(int i = 0; i < myTriangulation.getNumberOfCells(); i++){
  //     int triangleId = -1;
  //     int triangleNumber = myTriangulation.getCellTriangleNumber(i);
  //     printf("Cell #%d [", i);
  //
  //     for(int j = 0; j < 4; j++){
  //       int vertexId = -1;
  //       myTriangulation.getCellVertex(i, j, vertexId);
  //       printf(" %d", vertexId);
  //     }
  //     printf("]:\n");
  //
  //     for(int j = 0; j < triangleNumber; j++){
  //       myTriangulation.getCellTriangle(i, j, triangleId);
  //       printf("\ttriangle #%d [", triangleId);
  //
  //       for(int k = 0; k < 3; k++){
  //         int vertexId = -1;
  //         myTriangulation.getTriangleVertex(triangleId, k, vertexId);
  //         printf(" %d", vertexId);
  //       }
  //       printf("]\n");
  //     }
  //   }
  //
  //   myTriangulation.preprocessEdgeTriangles();
  //
  //
  //   for(int i = 0; i < myTriangulation.getNumberOfTriangles(); i++){
  //     vector<int> vertexIds(3);
  //     for(int j = 0; j < 3; j++){
  //       myTriangulation.getTriangleVertex(i, j, vertexIds[j]);
  //     }
  //     printf("Triangle #%d [%d, %d, %d]\n", i,
  //       vertexIds[0], vertexIds[1], vertexIds[2]);
  //   }
  //
  //   for(int i = 0; i < myTriangulation.getNumberOfEdges(); i++){
  //     int vertexId0, vertexId1;
  //     myTriangulation.getEdgeVertex(i, 0, vertexId0);
  //     myTriangulation.getEdgeVertex(i, 1, vertexId1);
  //     printf("Edge #%d [%d, %d]\n", i, vertexId0, vertexId1);
  //
  //     int triangleNumber = myTriangulation.getEdgeTriangleNumber(i);
  //     for(int j = 0; j < triangleNumber; j++){
  //       int triangleId;
  //       myTriangulation.getEdgeTriangle(i, j, triangleId);
  //       vector<int> vertexIds(3);
  //       for(int k = 0; k < 3; k++){
  //         myTriangulation.getTriangleVertex(triangleId, k, vertexIds[k]);
  //       }
  //
  //       printf("  triangle #%d [%d %d %d]\n",
  //         triangleId, vertexIds[0], vertexIds[1], vertexIds[2]);
  //     }
  //   }
  //
  //   myTriangulation.preprocessTriangleEdges();
  //
  //   for(int i = 0; i < myTriangulation.getNumberOfTriangles(); i++){
  //
  //     printf("triangle #%d\n\t", i);
  //
  //     for(int j = 0; j < 3; j++){
  //       int edgeId;
  //       myTriangulation.getTriangleEdge(i, j, edgeId);
  //       int vertexId0, vertexId1;
  //       myTriangulation.getEdgeVertex(edgeId, 0, vertexId0);
  //       myTriangulation.getEdgeVertex(edgeId, 1, vertexId1);
  //       printf("#%d [%d %d] ", edgeId, vertexId0, vertexId1);
  //     }
  //     printf("\n");
  //   }

  // end of test

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkScalarFieldNormalizer] Error: not enough input information."
         << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkScalarFieldNormalizer] Error: input has no point." << endl;
    return -1;
  }
#endif

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField = NULL;

  if(ScalarField.length()) {
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  } else {
    inputScalarField = input->GetPointData()->GetArray(0);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField) {
    cerr
      << "[ttkScalarFieldNormalizer] Error: input scalar field pointer is NULL."
      << endl;
    return -1;
  }
#endif

  if(outputScalarField_) {
    outputScalarField_->Delete();
    outputScalarField_ = NULL;
  }

  // allocate the memory for the output scalar field
  if(!outputScalarField_) {
    switch(inputScalarField->GetDataType()) {

      case VTK_CHAR:
        outputScalarField_ = vtkCharArray::New();
        break;

      case VTK_DOUBLE:
        outputScalarField_ = vtkDoubleArray::New();
        break;

      case VTK_FLOAT:
        outputScalarField_ = vtkFloatArray::New();
        break;

      case VTK_INT:
        outputScalarField_ = vtkIntArray::New();
        break;

      case VTK_ID_TYPE:
        outputScalarField_ = vtkIdTypeArray::New();
        break;

      default:
        stringstream msg;
        msg << "[ttkScalarFieldNormalizer] Unsupported data type :(" << endl;
        dMsg(cerr, msg.str(), fatalMsg);
        return -2;
    }
    outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
    outputScalarField_->SetName(inputScalarField->GetName());
  }

  // calling the executing package
  normalize(inputScalarField, outputScalarField_);

  // on the output, replace the field array by a pointer to its processed
  // version
  if(ScalarField.length()) {
    output->GetPointData()->RemoveArray(ScalarField.data());
  } else {
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);

  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkScalarFieldNormalizer::RequestData(vtkInformation *request,
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkScalarFieldNormalizer] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
