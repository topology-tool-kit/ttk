/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief command line program for mandatory critical point computation.

// include the local headers
#include <ttkMandatoryCriticalPoints.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

template <class vtkWriterClass>
int save(vtkDataSet *dataSet, stringstream &fileName) {

  vtkSmartPointer<vtkWriterClass> writer
    = vtkSmartPointer<vtkWriterClass>::New();

  writer->SetFileName(fileName.str().data());
  writer->SetInputData(dataSet);
  writer->Write();

  return 0;
}

int saveSaddles(vtkProgram<ttkMandatoryCriticalPoints> program,
                const int &portId,
                const int &saddleType) {

  vtkUnstructuredGrid *tree
    = (vtkUnstructuredGrid *)program.ttkObject_->GetOutput(portId);
  vtkDataArray *typeArray = tree->GetPointData()->GetArray("Type");
  vtkDataArray *componentArray = tree->GetPointData()->GetArray("ComponentId");

  for(int i = 0; i < (int)tree->GetNumberOfPoints(); i++) {
    double type = 0, componentId = 0;
    typeArray->GetTuple(i, &type);

    if(type == saddleType) {
      // this is a saddle
      componentArray->GetTuple(i, &componentId);

      string saddleTypeString = "join";
      if(saddleType != 1)
        saddleTypeString = "split";

      cout << "[main.cpp] Extracting " << saddleTypeString << " saddle #"
           << componentId << "..." << endl;
      program.ttkObject_->SetOutputJoinSaddleComponentId(componentId);
      program.ttkObject_->Update(saddleType);

      stringstream msg;
      msg << "output_" << saddleTypeString << "Saddle#" << componentId;

      if(program.ttkObject_->GetOutput(saddleType)) {
        if(program.ttkObject_->GetOutput(saddleType)->GetDataObjectType()
           == VTK_IMAGE_DATA) {
          msg << ".vti";
          save<vtkXMLImageDataWriter>(
            program.ttkObject_->GetOutput(saddleType), msg);
        }

        if(program.ttkObject_->GetOutput(saddleType)->GetDataObjectType()
           == VTK_POLY_DATA) {
          msg << ".vtp";
          save<vtkXMLPolyDataWriter>(
            program.ttkObject_->GetOutput(saddleType), msg);
        }

        if(program.ttkObject_->GetOutput(saddleType)->GetDataObjectType()
           == VTK_UNSTRUCTURED_GRID) {
          msg << ".vtu";
          save<vtkXMLUnstructuredGridWriter>(
            program.ttkObject_->GetOutput(saddleType), msg);
        }
      }
    }
  }

  return 0;
}

int main(int argc, char **argv) {

  // init editor
  vtkProgram<ttkMandatoryCriticalPoints> program;

  int lowerBoundId = 0, upperBoundId = 1;
  double normalizedSimplificationThreshold = 0;

  program.parser_.setArgument(
    "l", &lowerBoundId, "Point data field Id for the lower bound", true);
  program.parser_.setArgument(
    "u", &upperBoundId, "Point data field Id for the upper bound", true);
  program.parser_.setArgument(
    "s", &normalizedSimplificationThreshold,
    "Topological simplification threshold (normalized values)", true);

  int ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // setup the variables introduced above
  program.ttkObject_->SetLowerBoundField(lowerBoundId);
  program.ttkObject_->SetUpperBoundField(upperBoundId);
  program.ttkObject_->SetSimplificationThreshold(
    normalizedSimplificationThreshold);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  // HERE, exceptionnally, we will do more files (one per join saddle)
  saveSaddles(program, 4, 1);
  saveSaddles(program, 5, 2);

  return ret;
}
