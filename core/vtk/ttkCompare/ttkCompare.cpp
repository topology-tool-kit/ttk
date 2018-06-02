#include "ttkCompare.h"

#include <vtkCellData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCompare);

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
   if (meshOnly){
      output->ShallowCopy(input1);
   } else {
      // We override existing arrays
      output->DeepCopy(input1);
   }

   // clean data
   diffBFlag_  = 0;

   // computation
   computeMeshDiff(input1, input2, output);

   if (!meshOnly) {
      computeVertsDiff(input1, input2, output);
      computeCellsDiff(input1, input2, output);
   }

   addFlagFieldData(output);

   {
      stringstream msg;
      msg << "[ttkCompare] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}

// Private

void ttkCompare::computeMeshDiff(vtkDataSet *input1, vtkDataSet* input2, vtkDataSet *output)
{
   Triangulation *triangulation1 = ttkTriangulation::getTriangulation(input1);
   vtkSmartPointer<vtkUnsignedCharArray> diffVerts =
       createVTKArray<vtkUnsignedCharArray>("VerticesDiff", triangulation1->getNumberOfVertices());
   vtkSmartPointer<vtkUnsignedCharArray> diffCells =
       createVTKArray<vtkUnsignedCharArray>("CellsDiff", triangulation1->getNumberOfCells());

   unsigned char *vertArr = (unsigned char *)diffVerts->GetVoidPointer(0);
   unsigned char *cellArr = (unsigned char *)diffCells->GetVoidPointer(0);

   diffBFlag_ = diffBFlag_ | compare_.computeMeshDiff(vertArr, cellArr);

   output->GetPointData()->SetScalars(diffVerts);
   output->GetCellData()->AddArray(diffCells);
}

void ttkCompare::computeVertsDiff(vtkDataSet *input1, vtkDataSet* input2, vtkDataSet *output)
{
   const int nbVertsArray = output->GetPointData()->GetNumberOfArrays();
   // Compute the difference for each vertex array
   for (int va = 0; va < nbVertsArray; ++va) {
      vtkDataArray *inputScalarField1 = output->GetPointData()->GetArray(va);
      vtkDataArray *inputScalarField2 =
          input2->GetPointData()->GetArray(inputScalarField1->GetName());
      if (!inputScalarField1 || !inputScalarField2 ||
          inputScalarField1->GetNumberOfComponents() > 1 ||
          inputScalarField2->GetNumberOfComponents() > 1) {
         continue;
      }
      switch (inputScalarField1->GetDataType()) {
         vtkTemplateMacro({
            diffBFlag_ =
                diffBFlag_ | compare_.computeVertDiff<VTK_TT>(inputScalarField1->GetVoidPointer(0),
                                                              inputScalarField2->GetVoidPointer(0));
         });
      }
      output->GetPointData()->AddArray(inputScalarField1);
   }
}

void ttkCompare::computeCellsDiff(vtkDataSet *input1, vtkDataSet* input2, vtkDataSet *output)
{
   const int nbCellsArray = output->GetCellData()->GetNumberOfArrays();
   // Compute the difference for each cell array
   for (int va = 0; va < nbCellsArray; ++va) {
      vtkDataArray *inputScalarField1 = output->GetCellData()->GetArray(va);
      vtkDataArray *inputScalarField2 =
          input2->GetCellData()->GetArray(inputScalarField1->GetName());
      if (!inputScalarField1 || !inputScalarField2 ||
          inputScalarField1->GetNumberOfComponents() > 1 ||
          inputScalarField2->GetNumberOfComponents() > 1) {
         continue;
      }
      switch (inputScalarField1->GetDataType()) {
         vtkTemplateMacro({
            diffBFlag_ =
                diffBFlag_ | compare_.computeCellDiff<VTK_TT>(inputScalarField1->GetVoidPointer(0),
                                                              inputScalarField2->GetVoidPointer(0));
         });
      }
      output->GetCellData()->AddArray(inputScalarField1);
   }
}

void ttkCompare::addFlagFieldData(vtkDataSet *output)
{
   vtkSmartPointer<vtkUnsignedCharArray> bitFlagArray =
       createVTKArray<vtkUnsignedCharArray>("Compare", 1);
   if (debugLevel_ > 1) {
      std::cout << "[ttkCompare]: diff flag is " << diffBFlag_ << std::endl;
   }
   bitFlagArray->SetTuple1(0, diffBFlag_);
   output->GetFieldData()->AddArray(bitFlagArray);
}
