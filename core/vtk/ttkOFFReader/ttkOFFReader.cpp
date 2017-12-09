#include "ttkOFFReader.h"

#include "vtkCellType.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"

#include <iostream>
#include <sstream>

vtkStandardNewMacro(ttkOFFReader);

// Public
// {{{

void ttkOFFReader::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)")
     << "\n";
}

// }}}
// Protected
// {{{

ttkOFFReader::ttkOFFReader() {
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

int ttkOFFReader::RequestData(vtkInformation *request,
                           vtkInformationVector **inputVector,
                           vtkInformationVector *outputVector) {

  ifstream offFile(FileName, ios::in);

  if (!offFile) {
    cerr << "[ttkOFFReader]: Can't read file: " << FileName << endl;
    return -1;
  }

  std::string FileType;

  offFile >> FileType;
  if (FileType != "OFF") {
    cerr << "[ttkOFFReader]: Bad format for file: " << FileName << endl;
    return -2;
  }

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> mesh =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkDoubleArray> scalars =
      vtkSmartPointer<vtkDoubleArray>::New();

  int nbVerts, nbFaces, nbEdges, curLine;
  std::string line;

  offFile >> nbVerts >> nbFaces >> nbEdges;
  // go after line
  std::getline(offFile, line);

  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(nbVerts);
  scalars->SetName("Scalars");

  curLine = 0;
  while (curLine < nbVerts && std::getline(offFile, line)) {
    double x, y, z, w;
    std::stringstream ss(line);
    if (!(ss >> x >> y >> z >> w)) {
      // if the scalar is not on the file, we use an elevation on the x axis.
      w = x;
    }

    points->InsertNextPoint(x, y, z);
    scalars->SetTuple1(curLine, w);

    ++curLine;
  }

  mesh->SetPoints(points);
  mesh->GetPointData()->SetScalars(scalars);

  for (int i = 0; i < nbFaces; i++) {
    int nbCellVerts;
    vtkSmartPointer<vtkIdList> cellVerts = vtkSmartPointer<vtkIdList>::New();
    offFile >> nbCellVerts;
    for (int j = 0; j < nbCellVerts; j++) {
      int id;
      offFile >> id;
      cellVerts->InsertNextId(id);
    }
    switch (nbCellVerts) {
    case 2:
      mesh->InsertNextCell(VTK_LINE, cellVerts);
      break;
    case 3:
      mesh->InsertNextCell(VTK_TRIANGLE, cellVerts);
      break;
    case 4:
      mesh->InsertNextCell(VTK_TETRA, cellVerts);
      break;
    default:
      cerr << "[ttkOFFReader]: Unsupported cell type having: " << nbCellVerts
           << " vertices" << endl;
      return -3;
    }
  }

  cout << "[ttkOFFReader]: mesh have : " << mesh->GetNumberOfCells() << " cells"
       << endl;

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);

  // get the ouptut
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->ShallowCopy(mesh);

  return 1;
}

// }}}
