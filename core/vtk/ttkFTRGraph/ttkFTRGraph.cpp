#include "ttkFTRGraph.h"

vtkStandardNewMacro(ttkFTRGraph)

    int ttkFTRGraph::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs)
{
   Memory m;

   vtkDataSet *input  = inputs[0];
   vtkDataSet *output = outputs[0];

   Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

   if (!triangulation)
      return -1;

   triangulation->setWrapper(this);
   fTRGraph_.setupTriangulation(triangulation);
   fTRGraph_.setWrapper(this);

   // use a pointer-base copy for the input data -- to adapt if your wrapper does
   // not produce an output of the type of the input.
   output->ShallowCopy(input);

   // in the following, the target scalar field of the input is replaced in the
   // variable 'output' with the result of the computation.
   // if your wrapper produces an output of the same type of the input, you
   // should proceed in the same way.
   vtkDataArray *inputScalarField = NULL;

   if (ScalarField.length()) {
      inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
   } else {
      inputScalarField = input->GetPointData()->GetArray(0);
   }

   if (!inputScalarField)
      return -2;

   // allocate the memory for the output scalar field
   if (!outputScalarField_) {
      switch (inputScalarField->GetDataType()) {
         case VTK_CHAR:
            outputScalarField_ = vtkCharArray::New();
            break;

         case VTK_DOUBLE:
            outputScalarField_ = vtkDoubleArray::New();
            break;

         case VTK_FLOAT:
            outputScalarField_ = vtkFloatArray::New();
            break;

         case VTK_INT:
            outputScalarField_ = vtkIntArray::New();
            break;

            stringstream msg;
            msg << "[ttkFTRGraph] Unsupported data type :(" << endl;
            dMsg(cerr, msg.str(), fatalMsg);
      }
   }
   outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
   outputScalarField_->SetName(inputScalarField->GetName());

   // on the output, replace the field array by a pointer to its processed
   // version
   if (ScalarField.length()) {
      output->GetPointData()->RemoveArray(ScalarField.data());
   } else {
      output->GetPointData()->RemoveArray(0);
   }
   output->GetPointData()->AddArray(outputScalarField_);

   // calling the executing package
   switch (inputScalarField->GetDataType()) {
      vtkTemplateMacro({
         fTRGraph_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
         fTRGraph_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
         fTRGraph_.build<VTK_TT>(0);
      });
   }

   {
      stringstream msg;
      msg << "[ttkFTRGraph] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}
