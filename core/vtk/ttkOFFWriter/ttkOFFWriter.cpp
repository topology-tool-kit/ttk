#include <ttkOFFWriter.h>

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkOFFWriter);

// Public
// {{{

void ttkOFFWriter::PrintSelf(std::ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->Filename ? this->Filename : "(none)")
     << std::endl;
}

// }}}
// Protected
// {{{

ttkOFFWriter::ttkOFFWriter() {
  this->setDebugMsgPrefix("OFFWriter");
}

ttkOFFWriter::~ttkOFFWriter() {
}

int ttkOFFWriter::OpenFile() {

  std::ofstream f(Filename, ios::out);

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

  if(this->OpenFile() == -1) {
    this->printErr("Could not open file `" + std::string{FileName} + "' :(");
    return;
  }

  Stream << "OFF" << std::endl;
  Stream << dataSet->GetNumberOfPoints() << " " << dataSet->GetNumberOfCells()
         << " 0" << std::endl;

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

    Stream << std::endl;
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

    Stream << std::endl;
  }
}

// }}}
