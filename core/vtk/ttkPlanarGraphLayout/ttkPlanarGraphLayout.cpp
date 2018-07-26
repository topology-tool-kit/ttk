#include                  <ttkPlanarGraphLayout.h>

#include                  <unordered_set>
#include                  <boost/functional/hash.hpp>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPlanarGraphLayout)

int ttkPlanarGraphLayout::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

    Memory m;

    vtkUnstructuredGrid* input = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
    vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
    output->DeepCopy(input);

    vtkPoints* outputPointSet = (vtkPointSet::SafeDownCast(output))->GetPoints();


    planarGraphLayout_.setWrapper(this);

    Triangulation* triangulation = ttkTriangulation::getTriangulation(output);

    auto timeArray = input->GetFieldData()->GetAbstractArray("Time");
    auto branchArray = input->GetPointData()->GetAbstractArray("NodeBranch");
    auto branches = (size_t*) branchArray->GetVoidPointer(0);
    // auto sizeArray = input->GetPointData()->GetAbstractArray("Size");
    // auto size = (double*) sizeArray->GetVoidPointer(0);

    switch(outputPointSet->GetDataType()){

        vtkTemplateMacro(({
            planarGraphLayout_.execute<VTK_TT>(
                triangulation,
                branches,
                timeArray->GetSize(),
                (VTK_TT*) outputPointSet->GetVoidPointer(0)
            );
        }));
    }


    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Memory usage: " << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

  return 0;
}
