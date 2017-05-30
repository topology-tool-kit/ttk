#include                <ttkProgramBase.h>

int ttkProgramBase::execute(){

  if(!vtkWrapper_)
    return -1;
  
  for(int i = 0; i < (int) inputs_.size(); i++){
    vtkWrapper_->SetInputData(i, inputs_[i]);
  }
 
  vtkWrapper_->Update();
  vtkWrapper_->Modified();
  
  return 0;
}

template < class vtkReaderClass> 
  int ttkProgramBase::load(const string &fileName,
    vector<vtkSmartPointer<vtkReaderClass> > &readerList){
   
  readerList.resize(readerList.size() + 1);
  readerList.back() = vtkSmartPointer<vtkReaderClass>::New();
    
  readerList.back()->SetFileName(fileName.data());
  
  // handle debug messages
  {
    stringstream msg;
    msg << "[ttkProgramBase] Reading input data..." << endl;
    // choose where to display this message (cout, cerr, a file)
    // choose the priority of this message (1, nearly always displayed, 
    // higher values mean lower priorities)
    dMsg(cout, msg.str(), 1);
  }
  
  readerList.back()->Update();
  inputs_.push_back(readerList.back()->GetOutput());

  if(!inputs_.back())
    return -1;

  if(!inputs_.back()->GetNumberOfPoints())
    return -2;

  if(!inputs_.back()->GetNumberOfCells())
    return -3;

  {
    stringstream msg;
    msg << "[ttkProgramBase]   done! (read " 
      << inputs_.back()->GetNumberOfPoints()
      << " vertices, "
      << inputs_.back()->GetNumberOfCells() 
      << " cells)" << endl;
    dMsg(cout, msg.str(), Debug::infoMsg);
  }

  return 0;
}

int ttkProgramBase::load(const vector<string> &inputPaths){

  int ret = -1;

  for(int i = 0; i < (int) inputPaths.size(); i++){
    
    string extension = inputPaths[i].substr(
      inputPaths[i].find_last_of('.') + 1);
    
    if(extension == "vti"){
      ret = load<vtkXMLImageDataReader>(
        inputPaths[i], imageDataReaders_);
    }
    else if(extension == "vtp"){
      ret = load<vtkXMLPolyDataReader>(
        inputPaths[i], polyDataReaders_);
    }
    else if(extension == "vtu"){
      ret = load<vtkXMLUnstructuredGridReader>(
        inputPaths[i], unstructuredGridReaders_);
    }
    else{
      stringstream msg;
      msg << "[ttkProgramBase] Unkown input extension `" 
        << extension << "' :(" << endl;
      dMsg(cerr, msg.str(), Debug::fatalMsg);
      return -1;
    }

    if(ret)
      return ret;
  }

  return 0;
}

template <class vtkWriterClass>
  int ttkProgramBase::save(const int &outputPortId) const{
   
  if(!vtkWrapper_)
    return -1;
    
  string extension;
  
  if((vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType() 
    == VTK_IMAGE_DATA)
    ||(vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
      == TTK_IMAGE_DATA)){
    extension = "vti";
  }
  
  if((vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType() 
    == VTK_POLY_DATA)
    ||(vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
      == TTK_POLY_DATA)){
    extension = "vtp";
  }
  
  if((vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType() 
    == VTK_UNSTRUCTURED_GRID)
    ||(vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
      == TTK_UNSTRUCTURED_GRID)){
    extension = "vtu";
  }
  
  stringstream fileName;
  fileName << outputPath_
    << "_port#" << outputPortId << "." << extension;
  
  vtkSmartPointer<vtkWriterClass> writer = 
    vtkSmartPointer<vtkWriterClass>::New();
  writer->SetFileName(fileName.str().data());
  writer->SetInputData(vtkWrapper_->GetOutput(outputPortId));
  stringstream msg;
  msg << "[ttkProgramBase] Saving output file `" 
    << fileName.str() << "'..." << endl;
  dMsg(cout, msg.str(), Debug::infoMsg);
  
  writer->Write();
    
  return 0;
}

int ttkProgramBase::save() const{
 
  if(!vtkWrapper_)
    return -1;

  for(int i = 0; i < (int) vtkWrapper_->GetNumberOfOutputPorts(); i++){
    
    if(vtkWrapper_->GetOutput(i)){
  
      if((vtkWrapper_->GetOutput(i)->GetDataObjectType() == VTK_IMAGE_DATA)
        ||(vtkWrapper_->GetOutput(i)->GetDataObjectType() == TTK_IMAGE_DATA)){
        save<vtkXMLImageDataWriter>(i);
      }
      
      if((vtkWrapper_->GetOutput(i)->GetDataObjectType() == VTK_POLY_DATA)
        ||(vtkWrapper_->GetOutput(i)->GetDataObjectType() == TTK_POLY_DATA)){
        save<vtkXMLPolyDataWriter>(i);
      }
      
      if((vtkWrapper_->GetOutput(i)->GetDataObjectType() 
        == VTK_UNSTRUCTURED_GRID)
        ||(vtkWrapper_->GetOutput(i)->GetDataObjectType() ==
          TTK_UNSTRUCTURED_GRID)){
        save<vtkXMLUnstructuredGridWriter>(i);
      }
    }
  }
  
  return 0;
}
