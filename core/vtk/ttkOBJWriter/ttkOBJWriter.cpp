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

using namespace std;

vtkStandardNewMacro(ttkOBJWriter);

// Public
// {{{

void ttkOBJWriter::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->Filename ? this->Filename : "(none)")
     << endl;
}

// }}}
// Protected
// {{{

ttkOBJWriter::ttkOBJWriter() {
  Filename = NULL;
}

ttkOBJWriter::~ttkOBJWriter() {
  SetFilename(NULL);
}

int ttkOBJWriter::OpenFile() {

  ofstream f(Filename, ios::out);

  if(!f.fail()) {
    Stream = std::move(f);
  } else {
    return -1;
  }

  return 0;
}

void ttkOBJWriter::WriteData() {

  vtkDataSet *dataSet = vtkDataSet::SafeDownCast(this->GetInput());

  if(!dataSet)
    return;

  if(this->OpenFile()) {
    cerr << "[ttkOBJWriter] Could not open file `" << Filename << "' :("
         << endl;
    return;
  }

  double p[3];
  for(vtkIdType i = 0; i < dataSet->GetNumberOfPoints(); i++) {
    dataSet->GetPoint(i, p);
    Stream << "v " << p[0] << " " << p[1] << " " << p[2] << " " << endl;
  }

  for(vtkIdType i = 0; i < dataSet->GetNumberOfCells(); i++) {
    vtkCell *c = dataSet->GetCell(i);

    Stream << "f ";
    for(int j = 0; j < c->GetNumberOfPoints(); j++) {
      Stream << c->GetPointId(j) + 1 << " ";
    }

    Stream << endl;
  }
}

// }}}
