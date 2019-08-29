#include <ttkOFFWriter.h>

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

vtkStandardNewMacro(ttkOFFWriter);

// Public
// {{{

void ttkOFFWriter::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->Filename ? this->Filename : "(none)")
     << endl;
}

// }}}
// Protected
// {{{

ttkOFFWriter::ttkOFFWriter() {
  Filename = NULL;
}

ttkOFFWriter::~ttkOFFWriter() {
  SetFilename(NULL);
}

int ttkOFFWriter::OpenFile() {

  ofstream f(Filename, ios::out);

  if(!f.fail()) {
    Stream = std::move(f);
  } else {
    return -1;
  }

  return 0;
}

void ttkOFFWriter::WriteData() {

  vtkDataSet *dataSet = vtkDataSet::SafeDownCast(this->GetInput());

  if(!dataSet)
    return;

  if(this->OpenFile()) {
    cerr << "[ttkOFFWriter] Could not open file `" << FileName << "' :("
         << endl;
    return;
  }

  Stream << "OFF" << endl;
  Stream << dataSet->GetNumberOfPoints() << " " << dataSet->GetNumberOfCells()
         << " 0" << endl;

  double p[3];
  for(vtkIdType i = 0; i < dataSet->GetNumberOfPoints(); i++) {
    dataSet->GetPoint(i, p);
    Stream << p[0] << " " << p[1] << " " << p[2] << " ";

    // by default, store everything
    // use the field selector to select a subset
    for(int j = 0; j < dataSet->GetPointData()->GetNumberOfArrays(); j++) {
      vtkDataArray *array = dataSet->GetPointData()->GetArray(j);
      for(int k = 0; k < array->GetNumberOfComponents(); k++) {
        Stream << array->GetComponent(i, k) << " ";
      }
    }

    Stream << endl;
  }

  for(vtkIdType i = 0; i < dataSet->GetNumberOfCells(); i++) {
    vtkCell *c = dataSet->GetCell(i);

    Stream << c->GetNumberOfPoints() << " ";
    for(int j = 0; j < c->GetNumberOfPoints(); j++) {
      Stream << c->GetPointId(j) << " ";
    }

    // by default, store everything
    // use the field selector to select a subset
    for(int j = 0; j < dataSet->GetCellData()->GetNumberOfArrays(); j++) {
      vtkDataArray *array = dataSet->GetCellData()->GetArray(j);
      for(int k = 0; k < array->GetNumberOfComponents(); k++) {
        Stream << array->GetComponent(i, k) << " ";
      }
    }

    Stream << endl;
  }
}

// }}}
