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
   vtkSmartPointer<vtkCharArray> diffVerts =
       createVTKArray<vtkCharArray>("VerticesDiff", triangulation1->getNumberOfVertices());
   vtkSmartPointer<vtkCharArray> diffCells =
       createVTKArray<vtkCharArray>("CellsDiff", triangulation1->getNumberOfCells());
   compare_.setVertsArray(diffVerts->GetVoidPointer(0));
   compare_.setCellsArray(diffCells->GetVoidPointer(0));

   compare_.computeMeshDiff();

   output->GetPointData()->SetScalars(diffVerts);
   output->GetCellData()->AddArray(diffCells);

   if (!meshOnly) {
      std::cout << "Scalar difference not supported yet" << std::endl;
   }

   // TODO for each scalar (vert and cells)

   // in the following, the target scalar field of the input1 is replaced in the
   // variable 'output' with the result of the computation.
   // if your wrapper produces an output of the same type of the input1, you
   // should proceed in the same way.
   // vtkDataArray *inputScalarField = NULL;

   // if (ScalarField.length()) {
   //    inputScalarField = input1->GetPointData()->GetArray(ScalarField.data());
   // } else {
   //    inputScalarField = input1->GetPointData()->GetArray(0);
   // }

   // if (!inputScalarField)
   //    return -2;

   // allocate the memory for the output scalar field
   // if (!outputScalarField_) {
   //    switch (inputScalarField->GetDataType()) {
   //       case VTK_CHAR:
   //          outputScalarField_ = vtkCharArray::New();
   //          break;

   //       case VTK_DOUBLE:
   //          outputScalarField_ = vtkDoubleArray::New();
   //          break;

   //       case VTK_FLOAT:
   //          outputScalarField_ = vtkFloatArray::New();
   //          break;

   //       case VTK_INT:
   //          outputScalarField_ = vtkIntArray::New();
   //          break;

   //          stringstream msg;
   //          msg << "[ttkCompare] Unsupported data type :(" << endl;
   //          dMsg(cerr, msg.str(), fatalMsg);
   //    }
   // }
   // outputScalarField_->SetNumberOfTuples(input1->GetNumberOfPoints());
   // outputScalarField_->SetName(inputScalarField->GetName());

   // on the output, replace the field array by a pointer to its processed
   // version
   // if (ScalarField.length()) {
   //    output->GetPointData()->RemoveArray(ScalarField.data());
   // } else {
   //    output->GetPointData()->RemoveArray(0);
   // }
   // output->GetPointData()->AddArray(outputScalarField_);

   // calling the executing package
   // switch (inputScalarField->GetDataType()) {
   //    vtkTemplateMacro({
   //       compare_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
   //       compare_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
   //       compare_.execute<VTK_TT>(SomeIntegerArgument);
   //    });
   // }

   {
      stringstream msg;
      msg << "[ttkCompare] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}

