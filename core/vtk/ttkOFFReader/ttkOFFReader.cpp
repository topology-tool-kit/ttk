#include "ttkOFFReader.h"

#include "vtkCellType.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"

#include <iostream>
#include <sstream>

vtkStandardNewMacro(ttkOFFReader);

// Public
// {{{

void ttkOFFReader::PrintSelf(ostream &os, vtkIndent indent)
{
   this->Superclass::PrintSelf(os, indent);

   os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)") << endl;
}

// }}}
// Protected
// {{{

ttkOFFReader::ttkOFFReader()
{
   this->SetNumberOfInputPorts(0);
   this->SetNumberOfOutputPorts(1);
   this->FileName = NULL;
   this->nbFields_ = 0;
   this->mesh_     = vtkSmartPointer<vtkUnstructuredGrid>::New();
   this->points_   = vtkSmartPointer<vtkPoints>::New();
}

int ttkOFFReader::RequestData(vtkInformation *       request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector * outputVector)
{
   ifstream offFile(FileName, ios::in);

   if (!offFile) {
      cerr << "[ttkOFFReader] Can't read file: '" << FileName << "'" << endl;
      return -1;
   }

   std::string FileType;

   offFile >> FileType;
   if (FileType != "OFF") {
      cerr << "[ttkOFFReader] Bad format for file: '" 
        << FileName << "'" << endl;
      return -2;
   }

   int         nbVerts, nbFaces, curLine;
   std::string line;

   // init values
   offFile >> nbVerts >> nbFaces;
   std::getline(offFile, line);

   if (!nbVerts) {
      // empty file
      return 0;
   }

   curLine = 0;

   // count numbers of scalars
   std::getline(offFile, line);
   nbFields_ = countNbFields(line) - 3;

   // allocation
   scalars_.resize(nbFields_);
   for (int i = 0; i < nbFields_; i++) {
      scalars_[i] = vtkSmartPointer<vtkDoubleArray>::New();
      scalars_[i]->SetNumberOfComponents(1);
      scalars_[i]->SetNumberOfTuples(nbVerts);
      const std::string name = "ScalarField_" + std::to_string(i);
      scalars_[i]->SetName(name.c_str());
   }

   // process vertices
   while ((curLine = processLineVert(curLine, line)) <= nbVerts) {
      std::getline(offFile, line);
   }

   // add verts data to the mesh
   mesh_->SetPoints(points_);
   for (const auto &scalarArray : scalars_) {
      mesh_->GetPointData()->SetScalars(scalarArray);
   }

  if(line.length()){
    // don't miss the first simplex if any
    curLine = processLineCell(curLine, line);
  }

   // add cells to the mesh
   for (int i = 0; i < nbFaces; i++) {
      std::getline(offFile, line);
      curLine = processLineCell(curLine, line);
   }

#ifndef NDEBUG
   cout << "[ttkOFFReader] Read " << mesh_->GetNumberOfPoints() 
    << " vertice(s)" << endl;
   cout << "[ttkOFFReader] Read " << mesh_->GetNumberOfCells() 
    << " cell(s)" << endl;
#endif

   // get the info object
   vtkInformation *outInfo = outputVector->GetInformationObject(0);
   outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);

   // Set the output
   vtkUnstructuredGrid *output =
       vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   output->ShallowCopy(mesh_);

   return 1;
}

int ttkOFFReader::countNbFields(std::string line)
{
   std::istringstream ss(line);
   double            buffer;
   int               nbLine = -1;
   while (!ss.fail()) {
      ss >> buffer;
      ++nbLine;
   }
   return nbLine;
}

int ttkOFFReader::processLineVert(int curLine, std::string &line)
{
   double             x, y, z;
   std::istringstream ss(line);

   // Coords
   ss >> x >> y >> z;

   points_->InsertNextPoint(x, y, z);

   // Scalars
   for (int i = 0; i < nbFields_; i++) {
      double scalar;
      ss >> scalar;
      scalars_[i]->SetTuple1(curLine, scalar);
   }

   return ++curLine;
}

int ttkOFFReader::processLineCell(int curLine, std::string &line)
{
   int                        nbCellVerts;
   vtkSmartPointer<vtkIdList> cellVerts = vtkSmartPointer<vtkIdList>::New();
   std::istringstream ss(line);
   ss >> nbCellVerts;
   for (int j = 0; j < nbCellVerts; j++) {
      int id;
      ss >> id;
      cellVerts->InsertNextId(id);
   }
   switch (nbCellVerts) {
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
         cerr << "[ttkOFFReader] Unsupported cell type having " << nbCellVerts << " vertices"
              << endl;
         return -3;
   }

   return ++curLine;
}

// }}}
