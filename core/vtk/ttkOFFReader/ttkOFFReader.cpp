#include <ttkOFFReader.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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
  this->setDebugMsgPrefix("OFFReader");
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

int countVertsData(const std::string &line) {
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

int countCellsData(const std::string &line) {
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

int processLineVert(vtkIdType curLine,
                    std::string &line,
                    vtkPoints *points,
                    std::vector<vtkNew<vtkDoubleArray>> &vertScalars,
                    const vtkIdType nbVertsData) {

  double x, y, z;
  std::istringstream ss(line);

  // Coords
  ss >> x >> y >> z;

  points->InsertNextPoint(x, y, z);

  // Scalars data
  for(vtkIdType i = 0; i < nbVertsData; i++) {
    double scalar;
    ss >> scalar;
    vertScalars[i]->SetTuple1(curLine, scalar);
  }

  return ++curLine;
}

int processLineCell(vtkIdType curLine,
                    std::string &line,
                    vtkUnstructuredGrid *mesh,
                    std::vector<vtkNew<vtkDoubleArray>> &cellScalars,
                    const vtkIdType nbVerts,
                    const vtkIdType nbCellsData,
                    const ttk::Debug &dbg) {
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
      mesh->InsertNextCell(VTK_LINE, cellVerts);
      break;
    case 3:
      mesh->InsertNextCell(VTK_TRIANGLE, cellVerts);
      break;
    case 4:
      mesh->InsertNextCell(VTK_TETRA, cellVerts);
      break;
    default:
      dbg.printErr("Unsupported cell type having " + std::to_string(nbCellVerts)
                   + " vertices");
      return -3;
  }

  // Scalars data
  for(vtkIdType i = 0; i < nbCellsData; i++) {
    double scalar;
    ss >> scalar;
    // Currline is after having read all the vertices
    cellScalars[i]->SetTuple1(curLine - nbVerts, scalar);
  }

  return ++curLine;
}

int ttkOFFReader::RequestData(vtkInformation *ttkNotUsed(request),
                              vtkInformationVector **ttkNotUsed(inputVector),
                              vtkInformationVector *outputVector) {
  std::ifstream offFile(FileName, ios::in);

  if(!offFile) {
    this->printErr("Can't read file: '" + std::string{FileName} + "'");
    return -1;
  }

  std::string FileType;

  offFile >> FileType;
  if(FileType != "OFF") {
    this->printErr("Bad format for file: '" + std::string{FileName} + "'");
    return -2;
  }

  vtkIdType curLine;
  std::string line;
  vtkIdType nbVerts{}, nbCells{}, nbVertsData{}, nbCellsData{};

  // init values
  offFile >> nbVerts >> nbCells;
  std::getline(offFile, line);

  if(nbVerts == 0) {
    // empty file
    return 0;
  }

  curLine = 0;

  // count numbers of vertices scalars
  std::getline(offFile, line);
  nbVertsData = countVertsData(line);

  // allocation verts
  std::vector<vtkNew<vtkDoubleArray>> vertScalars{};
  vertScalars.resize(nbVertsData);
  for(vtkIdType i = 0; i < nbVertsData; i++) {
    vertScalars[i]->SetNumberOfComponents(1);
    vertScalars[i]->SetNumberOfTuples(nbVerts);
    const std::string name = "VertScalarField_" + std::to_string(i);
    vertScalars[i]->SetName(name.c_str());
  }

  vtkNew<vtkPoints> points{};

  // read vertices
  while(
    (curLine = processLineVert(curLine, line, points, vertScalars, nbVertsData))
    < nbVerts) {
    std::getline(offFile, line);
  }

  // add verts data to the mesh
  vtkNew<vtkUnstructuredGrid> mesh{};
  mesh->SetPoints(points);
  for(const auto &scalarArray : vertScalars) {
    mesh->GetPointData()->AddArray(scalarArray);
  }

  // count numbers of cells scalars
  std::getline(offFile, line);
  nbCellsData = countCellsData(line);

  // allocation cells
  std::vector<vtkNew<vtkDoubleArray>> cellScalars{};
  cellScalars.resize(nbCellsData);
  for(vtkIdType i = 0; i < nbCellsData; i++) {
    cellScalars[i]->SetNumberOfComponents(1);
    cellScalars[i]->SetNumberOfTuples(nbCells);
    const std::string name = "CellScalarField_" + std::to_string(i);
    cellScalars[i]->SetName(name.c_str());
  }

  // read cells
  if(nbCells != 0) {
    while((curLine = processLineCell(
             curLine, line, mesh, cellScalars, nbVerts, nbCellsData, *this))
          < nbVerts + nbCells) {
      std::getline(offFile, line);
    }
  }

  // add cell data to the mesh
  for(const auto &scalarArray : cellScalars) {
    mesh->GetCellData()->AddArray(scalarArray);
  }

#ifndef NDEBUG
  this->printMsg("Read " + std::to_string(mesh->GetNumberOfPoints())
                 + " vertice(s)");
  this->printMsg("Read " + std::to_string(mesh->GetNumberOfCells())
                 + " cell(s)");
#endif

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);

  // Set the output
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->ShallowCopy(mesh);

  return 1;
}

// }}}
