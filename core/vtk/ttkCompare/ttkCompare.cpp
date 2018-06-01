#include "ttkCompare.h"

#include <vtkCellData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCompare)

    int ttkCompare::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs)
{
   Memory m;

   vtkDataSet *   input1         = inputs[0];
   vtkDataSet *   input2         = inputs[1];
   Triangulation *triangulation1 = ttkTriangulation::getTriangulation(input1);
   Triangulation *triangulation2 = ttkTriangulation::getTriangulation(input2);

   if (!triangulation1 || !triangulation2)
      return -1;

   triangulation1->setWrapper(this);
   compare_.setupTriangulations(triangulation1, triangulation2);
   compare_.setWrapper(this);

   vtkDataSet *output = outputs[0];
   output->ShallowCopy(input1);

   // Store difference for mesh
   vtkSmartPointer<vtkUnsignedCharArray> diffVerts =
       createVTKArray<vtkUnsignedCharArray>("VerticesDiff", triangulation1->getNumberOfVertices());
   vtkSmartPointer<vtkUnsignedCharArray> diffCells =
       createVTKArray<vtkUnsignedCharArray>("CellsDiff", triangulation1->getNumberOfCells());
   compare_.setVertsArray(diffVerts->GetVoidPointer(0));
   compare_.setCellsArray(diffCells->GetVoidPointer(0));

   diffReturn_ = compare_.computeMeshDiff();

   output->GetPointData()->SetScalars(diffVerts);
   output->GetCellData()->AddArray(diffCells);

   if (!meshOnly) {
      const int nbVertsArray = input1->GetPointData()->GetNumberOfArrays();
      const int nbCellsArray = input1->GetCellData()->GetNumberOfArrays();

      // Remove other scalars to replace these arrays with their compared counterpart
      for(int va = 0; va < nbVertsArray; ++va) {
         output->GetPointData()->RemoveArray(0);
      }
      for(int ca = 0; ca < nbCellsArray; ++ca) {
         output->GetCellData()->RemoveArray(0);
      }

      // Compute the difference for each array
      for(int va = 0; va < nbVertsArray; ++va) {
         vtkDataArray *inputScalarField1 = nullptr;
         vtkDataArray *inputScalarField2 = nullptr;
         inputScalarField1               = input1->GetPointData()->GetArray(va);
         inputScalarField2 = input2->GetPointData()->GetArray(inputScalarField1->GetName());
         if (!inputScalarField1 || !inputScalarField2 ||
             inputScalarField1->GetNumberOfComponents() > 1 ||
             inputScalarField2->GetNumberOfComponents() > 1) {
            continue;
         }
         switch (inputScalarField1->GetDataType()) {
            vtkTemplateMacro({
               compare_.computeVertDiff<VTK_TT>(inputScalarField1->GetVoidPointer(0),
                                                  inputScalarField2->GetVoidPointer(0));
            });
         }
         output->GetPointData()->AddArray(inputScalarField1);
      }

      // Compute the difference for each array
      for (int va = 0; va < nbCellsArray; ++va) {
         vtkDataArray *inputScalarField1 = nullptr;
         vtkDataArray *inputScalarField2 = nullptr;
         inputScalarField1               = input1->GetCellData()->GetArray(va);
         inputScalarField2 = input2->GetCellData()->GetArray(inputScalarField1->GetName());
         if (!inputScalarField1 || !inputScalarField2 ||
             inputScalarField1->GetNumberOfComponents() > 1 ||
             inputScalarField2->GetNumberOfComponents() > 1) {
            continue;
         }
         switch (inputScalarField1->GetDataType()) {
            vtkTemplateMacro({
               compare_.computeCellDiff<VTK_TT>(inputScalarField1->GetVoidPointer(0),
                                                inputScalarField2->GetVoidPointer(0));
            });
         }
         output->GetCellData()->AddArray(inputScalarField1);
      }
   }

   {
      stringstream msg;
      msg << "[ttkCompare] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}

