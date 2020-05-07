#include <ttkScalarFieldNormalizer.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkScalarFieldNormalizer);

ttkScalarFieldNormalizer::ttkScalarFieldNormalizer() {
  this->setDebugMsgPrefix("ScalarFieldNormalizer");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkScalarFieldNormalizer::~ttkScalarFieldNormalizer() {
}

int ttkScalarFieldNormalizer::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkScalarFieldNormalizer::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkScalarFieldNormalizer::normalize(vtkDataArray *input,
                                        vtkDataArray *output) const {
  if(!output)
    return -1;
  if(!input)
    return -2;

  double min = 0, max = 0;
  for(ttk::SimplexId i = 0; i < input->GetNumberOfTuples(); i++) {

    double value = input->GetTuple1(i);

    if((!i) || (value < min)) {
      min = value;
    }
    if((!i) || (value > max)) {
      max = value;
    }
  }

  for(ttk::SimplexId i = 0; i < input->GetNumberOfTuples(); i++) {
    double value = input->GetTuple1(i);

    value = (value - min) / (max - min) + pow10(-FLT_DIG);

    output->SetTuple1(i, value);
  }

  return 0;
}

int ttkScalarFieldNormalizer::RequestData(vtkInformation *request,
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {

  ttk::Timer t;

  this->printMsg("Normalizing Scalar Field", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  // test
  //   ttkTriangulation myTriangulation;
  //   myTriangulation.setInputData(input);
  //
  //   myTriangulation.preconditionCellTriangles();
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
  //   myTriangulation.preconditionEdgeTriangles();
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
  //   myTriangulation.preconditionTriangleEdges();
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
    this->printErr("Not enough input information.");
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
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
  auto inputScalarField = this->GetInputArrayToProcess(0, inputVector);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField) {
    this->printErr("Input scalar field pointer is NULL.");
    return -1;
  }
#endif

  auto outputScalarField
    = vtkSmartPointer<vtkDataArray>::Take(inputScalarField->NewInstance());
  outputScalarField->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField->SetName(inputScalarField->GetName());

  // calling the executing package
  normalize(inputScalarField, outputScalarField);

  auto inputArrayAssociation = this->GetInputArrayAssociation(0, inputVector);
  if(inputArrayAssociation == 0) {
    output->GetPointData()->AddArray(outputScalarField);
  } else if(inputArrayAssociation == 1) {
    output->GetCellData()->AddArray(outputScalarField);
  } else {
    output->GetFieldData()->AddArray(outputScalarField);
  }

  this->printMsg(
    "Normalizing Scalar Field", 1, t.getElapsedTime(), this->threadNumber_);

  return 1;
}