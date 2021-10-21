#include <Triangulation.h>
#include <ttkTriangulationReader.h>

#include <vtkDataObject.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <ttkUtils.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTriangulationReader);

ttkTriangulationReader::ttkTriangulationReader() {
  this->setDebugMsgPrefix("TriangulationReader");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTriangulationReader::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkTriangulationReader::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkTriangulationReader::validateFilePath() {
  if(this->TriangulationFilePath.length() < 8
     || this->TriangulationFilePath
            .substr(this->TriangulationFilePath.length() - 4, 4)
            .compare(".tpt")
          != 0) {
    this->printErr(
      "TTK Preconditioned Triangulation file has to end with '.tpt'.");
    return 0;
  }

  return 1;
}

int ttkTriangulationReader::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  ttk::Timer timer;

  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector, 0);

  if(input == nullptr || output == nullptr) {
    this->printErr("Invalid input dataset");
    return 0;
  }

  output->ShallowCopy(input);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Invalid input triangulation");
    return 0;
  }

  auto explTri
    = static_cast<ttk::ExplicitTriangulation *>(triangulation->getData());

  if(!this->validateFilePath()) {
    return 0;
  }
  std::ifstream in(this->TriangulationFilePath);
  explTri->readFromFile(in);

  this->printMsg("Restored triangulation from " + this->TriangulationFilePath,
                 1.0, timer.getElapsedTime(), 1);

  return 1;
}
