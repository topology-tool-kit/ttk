#include <ttkOBJWriter.h>

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

#include <iostream>
#include <sstream>

vtkStandardNewMacro(ttkOBJWriter);

// Public
// {{{

void ttkOBJWriter::PrintSelf(ostream &os, vtkIndent indent){
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " 
    << (this->FileName ? this->FileName : "(none)") << endl;
}

// }}}
// Protected
// {{{

ttkOBJWriter::ttkOBJWriter(){
  FileName = NULL;
  Stream = NULL;
}

ttkOBJWriter::~ttkOBJWriter(){
  SetFileName(NULL);
  if(Stream)
    delete Stream;
}

int ttkOBJWriter::OpenFile(){

  ofstream *f = new ofstream(FileName, ios::out);
  
  if(!f->fail()){
    Stream = f;
  }
  else{
    delete f;
    return -1;
  }
  
  return 0;
}

void ttkOBJWriter::WriteData(){
  
  vtkDataSet *dataSet = 
    vtkDataSet::SafeDownCast(this->GetInput());
    
  if(!dataSet)
    return;
  
  if(this->OpenFile()){
    cerr << "[ttkOBJWriter] Could not open file `"
      << FileName << "' :(" << endl; 
    return;
  }
  
  double p[3];
  for(int i = 0; i < dataSet->GetNumberOfPoints(); i++){
    dataSet->GetPoint(i, p);
    (*Stream) << "v " <<  p[0] << " " << p[1] << " " << p[2] << " " << endl;
  }
  
  for(int i = 0; i < dataSet->GetNumberOfCells(); i++){
    vtkCell *c = dataSet->GetCell(i);
   
    (*Stream) << "f ";
    for(int j = 0; j < c->GetNumberOfPoints(); j++){
      (*Stream) << c->GetPointId(j) + 1 << " ";
    }
    
    (*Stream) << endl;
  }
}


// }}}
