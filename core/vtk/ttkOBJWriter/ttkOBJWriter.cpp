#include <ttkOBJWriter.h>

#include <vtkCell.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkOBJWriter);

// Public
// {{{

void ttkOBJWriter::PrintSelf(std::ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->Filename ? this->Filename : "(none)")
     << std::endl;
}

// }}}
// Protected
// {{{

ttkOBJWriter::ttkOBJWriter() {
  this->setDebugMsgPrefix("OBJWriter");
}

ttkOBJWriter::~ttkOBJWriter() {
}

int ttkOBJWriter::OpenFile() {

  std::ofstream f(Filename, ios::out);

  if(!f.fail()) {
    Stream = std::move(f);
  } else {
    return -1;
  }

  return 0;
}

void ttkOBJWriter::WriteData() {

  auto dataSet = vtkDataSet::SafeDownCast(this->GetInput());

  if(dataSet == nullptr)
    return;

  if(this->OpenFile() == -1) {
    this->printErr("Could not open file `" + std::string{Filename} + "' :(");
    return;
  }

  double p[3];
  for(vtkIdType i = 0; i < dataSet->GetNumberOfPoints(); i++) {
    dataSet->GetPoint(i, p);
    Stream << "v " << p[0] << " " << p[1] << " " << p[2] << " " << std::endl;
  }

  for(vtkIdType i = 0; i < dataSet->GetNumberOfCells(); i++) {
    vtkCell *c = dataSet->GetCell(i);

    Stream << "f ";
    for(int j = 0; j < c->GetNumberOfPoints(); j++) {
      Stream << c->GetPointId(j) + 1 << " ";
    }

    Stream << std::endl;
  }
}

// }}}
