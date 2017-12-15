#include <ttkOFFWriter.h>

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

vtkStandardNewMacro(ttkOFFWriter);

// Public
// {{{

void ttkOFFWriter::PrintSelf(ostream &os, vtkIndent indent){
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " 
    << (this->FileName ? this->FileName : "(none)") << endl;
}

// }}}
// Protected
// {{{

ttkOFFWriter::ttkOFFWriter(){
  FileName = NULL;
  Stream = NULL;
}

ttkOFFWriter::~ttkOFFWriter(){
  SetFileName(NULL);
  if(Stream)
    delete Stream;
}

void ttkOFFWriter::WriteData(){
  
  vtkDataSet *grid = 
    vtkDataSet::SafeDownCast(this->GetInput());
    
  if(!grid)
    return;
  
  cout << "Writing things!!!!" << endl;
}


// }}}
