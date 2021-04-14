#include <Triangulation.h>
#include <ttkTriangulationWriter.h>

#include <array>

#include <vtkDataArray.h>
#include <vtkExecutive.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTriangulationWriter);

ttkTriangulationWriter::ttkTriangulationWriter() {
  this->SetNumberOfInputPorts(1);
  this->setDebugMsgPrefix("TriangulationWriter");
}

int ttkTriangulationWriter::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTriangulationWriter::OpenFile() {

  std::ofstream f(Filename, std::ios::out | std::ios::binary);

  if(!f.fail()) {
    Stream = std::move(f);
  } else {
    return -1;
  }

  return 0;
}

template <typename T>
void writeBin(std::ofstream &stream, const T var) {
  stream.write(reinterpret_cast<const char *>(&var), sizeof(var));
}

int ttkTriangulationWriter::Write() {

  ttk::Timer tm{};
  this->setDebugLevel(3);

  const auto dataset = vtkDataSet::SafeDownCast(this->GetInput());
  if(dataset == nullptr) {
    this->printErr("Invalid dataset");
    return 1;
  }
  if(!dataset->IsA("vtkUnstructuredGrid")) {
    this->printErr("This filter only works on Explicit Triangulations");
    return 1;
  }
  const auto triangulation = ttkAlgorithm::GetTriangulation(dataset);
  if(triangulation == nullptr) {
    this->printErr("Invalid triangulation");
    return 1;
  }

  const auto explTri
    = static_cast<ttk::ExplicitTriangulation *>(triangulation->getData());

  if(this->OpenFile() == -1) {
    this->printErr("Could not open file `" + std::string{Filename} + "' :(");
    return 1;
  }

  explTri->writeToFile(this->Stream);
  this->Stream.flush();

  this->printMsg("Wrote triangulation to " + std::string{this->Filename}, 1.0,
                 tm.getElapsedTime(), 1);

  return 0;
}

vtkDataObject *ttkTriangulationWriter::GetInput() {
  // copied from ParaView's vtkWriter::GetInput()
  if(this->GetNumberOfInputConnections(0) < 1) {
    return nullptr;
  }
  return this->GetExecutive()->GetInputData(0, 0);
}

void ttkTriangulationWriter::SetInputData(vtkDataObject *input) {
  // copied from ParaView's vtkWriter::SetInputData()
  this->SetInputDataInternal(0, input);
}
