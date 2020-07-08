#include "ttkOFFReader.h"

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"
#include <vtkDataSetReader.h>

#include <iostream>
#include <sstream>

vtkStandardNewMacro(ttkOFFReader);

// Public
// {{{

void ttkOFFReader::PrintSelf(std::ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)")
     << std::endl;
}

// }}}
// Protected
// {{{

ttkOFFReader::ttkOFFReader() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

int ttkOFFReader::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
  std::ifstream offFile(FileName, ios::in);

  if(!offFile) {
    std::cerr << "[ttkOFFReader] Can't read file: '" << FileName << "'"
              << std::endl;
    return -1;
  }

  std::string FileType;

  offFile >> FileType;
  if(FileType != "OFF") {
    std::cerr << "[ttkOFFReader] Bad format for file: '" << FileName << "'"
              << std::endl;
    return -2;
  }

  vtkIdType curLine;
  std::string line;

  // init values
  offFile >> nbVerts_ >> nbCells_;
  std::getline(offFile, line);

  if(!nbVerts_) {
    // empty file
    return 0;
  }

  curLine = 0;

  // count numbers of vertices scalars
  std::getline(offFile, line);
  nbVertsData_ = countVertsData(line);

  // allocation verts
  vertScalars_.resize(nbVertsData_);
  for(vtkIdType i = 0; i < nbVertsData_; i++) {
    vertScalars_[i]->SetNumberOfComponents(1);
    vertScalars_[i]->SetNumberOfTuples(nbVerts_);
    const std::string name = "VertScalarField_" + std::to_string(i);
    vertScalars_[i]->SetName(name.c_str());
  }

  // read vertices
  while((curLine = processLineVert(curLine, line)) < nbVerts_) {
    std::getline(offFile, line);
  }

  // add verts data to the mesh
  mesh_->SetPoints(points_);
  for(const auto &scalarArray : vertScalars_) {
    mesh_->GetPointData()->AddArray(scalarArray);
  }

  // count numbers of cells scalars
  std::getline(offFile, line);
  nbCellsData_ = countCellsData(line);

  // allocation cells
  cellScalars_.resize(nbCellsData_);
  for(vtkIdType i = 0; i < nbCellsData_; i++) {
    cellScalars_[i]->SetNumberOfComponents(1);
    cellScalars_[i]->SetNumberOfTuples(nbCells_);
    const std::string name = "CellScalarField_" + std::to_string(i);
    cellScalars_[i]->SetName(name.c_str());
  }

  // read cells
  if(nbCells_) {
    while((curLine = processLineCell(curLine, line)) < nbVerts_ + nbCells_) {
      std::getline(offFile, line);
    }
  }

  // add cell data to the mesh
  for(const auto &scalarArray : cellScalars_) {
    mesh_->GetCellData()->AddArray(scalarArray);
  }

#ifndef NDEBUG
  std::cout << "[ttkOFFReader] Read " << mesh_->GetNumberOfPoints()
            << " vertice(s)" << std::endl;
  std::cout << "[ttkOFFReader] Read " << mesh_->GetNumberOfCells() << " cell(s)"
            << std::endl;
#endif

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);

  // Set the output
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->ShallowCopy(mesh_);

  return 1;
}

int ttkOFFReader::countVertsData(std::string line) {
  std::istringstream ss(line);
  double buffer;
  int nbFields = -1;
  while(!ss.fail()) {
    ss >> buffer;
    ++nbFields;
  }
  // remove the 3 coordinates
  return nbFields - 3;
}

int ttkOFFReader::countCellsData(std::string line) {
  std::istringstream ss(line);
  double buffer;
  int nbFields = -1;
  int sizeCell;
  // This line start by the size of the cell
  ss >> sizeCell;
  // we read all the fields
  while(!ss.fail()) {
    ss >> buffer;
    ++nbFields;
  }
  // remove the vertices of the cell
  return nbFields - sizeCell;
}

int ttkOFFReader::processLineVert(vtkIdType curLine, std::string &line) {

  double x, y, z;
  std::istringstream ss(line);

  // Coords
  ss >> x >> y >> z;

  points_->InsertNextPoint(x, y, z);

  // Scalars data
  for(vtkIdType i = 0; i < nbVertsData_; i++) {
    double scalar;
    ss >> scalar;
    vertScalars_[i]->SetTuple1(curLine, scalar);
  }

  return ++curLine;
}

int ttkOFFReader::processLineCell(vtkIdType curLine, std::string &line) {
  int nbCellVerts;
  vtkNew<vtkIdList> cellVerts{};
  std::istringstream ss(line);
  ss >> nbCellVerts;
  for(int j = 0; j < nbCellVerts; j++) {
    vtkIdType id;
    ss >> id;
    cellVerts->InsertNextId(id);
  }
  switch(nbCellVerts) {
    case 2:
      mesh_->InsertNextCell(VTK_LINE, cellVerts);
      break;
    case 3:
      mesh_->InsertNextCell(VTK_TRIANGLE, cellVerts);
      break;
    case 4:
      mesh_->InsertNextCell(VTK_TETRA, cellVerts);
      break;
    default:
      std::cerr << "[ttkOFFReader] Unsupported cell type having " << nbCellVerts
                << " vertices" << std::endl;
      return -3;
  }

  // Scalars data
  for(vtkIdType i = 0; i < nbCellsData_; i++) {
    double scalar;
    ss >> scalar;
    // Currline is after having read all the vertices
    cellScalars_[i]->SetTuple1(curLine - nbVerts_, scalar);
  }

  return ++curLine;
}

// }}}
