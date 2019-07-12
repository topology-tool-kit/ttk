#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int ttkProgramBase::execute() {

  if(!vtkWrapper_)
    return -1;

  for(int i = 0; i < (int)inputs_.size(); i++) {
    vtkWrapper_->SetInputData(i, inputs_[i]);
  }

  vtkWrapper_->Update();
  vtkWrapper_->Modified();

  return 0;
}

int ttkProgramBase::load(const vector<string> &inputPaths) {

  int ret = -1;

  for(int i = 0; i < (int)inputPaths.size(); i++) {

    string extension
      = inputPaths[i].substr(inputPaths[i].find_last_of('.') + 1);

    if(extension == "vti") {
      ret = load<vtkXMLImageDataReader>(inputPaths[i], imageDataReaders_);
    } else if(extension == "vtp") {
      ret = load<vtkXMLPolyDataReader>(inputPaths[i], polyDataReaders_);
    } else if(extension == "vtu") {
      ret = load<vtkXMLUnstructuredGridReader>(
        inputPaths[i], unstructuredGridReaders_);
    } else {
      stringstream msg;
      msg << "[ttkProgramBase] Unkown input extension `" << extension << "' :("
          << endl;
      dMsg(cerr, msg.str(), Debug::fatalMsg);
      return -1;
    }

    if(ret)
      return ret;
  }

  return 0;
}

int ttkProgramBase::save() const {

  if(!vtkWrapper_)
    return -1;

  for(int i = 0; i < vtkWrapper_->GetNumberOfOutputPorts(); i++) {

    if(vtkWrapper_->GetOutput(i)) {

      if((vtkWrapper_->GetOutput(i)->GetDataObjectType() == VTK_IMAGE_DATA)) {
        //         ||(vtkWrapper_->GetOutput(i)->GetDataObjectType() ==
        //         TTK_IMAGE_DATA)){
        save<vtkXMLImageDataWriter>(i);
      }

      if((vtkWrapper_->GetOutput(i)->GetDataObjectType() == VTK_POLY_DATA)) {
        //         ||(vtkWrapper_->GetOutput(i)->GetDataObjectType() ==
        //         TTK_POLY_DATA)){
        save<vtkXMLPolyDataWriter>(i);
      }

      if((vtkWrapper_->GetOutput(i)->GetDataObjectType()
          == VTK_UNSTRUCTURED_GRID)) {
        //         ||(vtkWrapper_->GetOutput(i)->GetDataObjectType() ==
        //           TTK_UNSTRUCTURED_GRID)){
        save<vtkXMLUnstructuredGridWriter>(i);
      }
    }
  }

  return 0;
}
