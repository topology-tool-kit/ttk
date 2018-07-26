#include                  <ttkDepthImageBasedGeometryApproximation.h>

using namespace std;
using namespace ttk;

#include <vtkXMLImageDataWriter.h>

vtkStandardNewMacro(ttkDepthImageBasedGeometryApproximation)

int ttkDepthImageBasedGeometryApproximation::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

    Memory m;

    vtkUnstructuredGrid* outputMesh =  vtkUnstructuredGrid::SafeDownCast(outputs[0]);
    vtkImageData* inputImage =  vtkImageData::SafeDownCast(inputs[0]);

    depthImageBasedGeometryApproximation_.setWrapper(this);

    auto depthValues = inputImage->GetPointData()->GetAbstractArray("DepthValues");

    auto camHeight = inputImage->GetFieldData()->GetAbstractArray("CamHeight");
    auto camPosition = inputImage->GetFieldData()->GetAbstractArray("CamPosition");
    auto camDirection = inputImage->GetFieldData()->GetAbstractArray("CamDirection");
    auto camUp = inputImage->GetFieldData()->GetAbstractArray("CamUp");
    auto camNearFar = inputImage->GetFieldData()->GetAbstractArray("CamNearFar");
    auto camRes = inputImage->GetFieldData()->GetAbstractArray("CamRes");

    if(depthValues==nullptr || camHeight==nullptr || camPosition==nullptr || camDirection==nullptr || camUp==nullptr || camNearFar==nullptr || camRes==nullptr){
        stringstream msg;
        msg << "[ttkDepthImageBasedGeometryApproximation] ERROR: Input depth image does not have one or more of the required fields (see SpecX - Appendix)" << endl;
        dMsg(cout, msg.str(), memoryMsg);
        return -1;
    }

    // depthImageBasedGeometryApproximation_.setInputDataPointer(depthValues->GetVoidPointer(0));

    vector<tuple<float,float,float>> vertices;
    vector<tuple<int,int,int>> triangles;
    vector<float> triangleDistortions;

    switch(depthValues->GetDataType()){
        vtkTemplateMacro({
            depthImageBasedGeometryApproximation_.setInputDataPointer( depthValues->GetVoidPointer(0) );

            depthImageBasedGeometryApproximation_.execute<VTK_TT>(
                // Arguments
                (double*) camPosition->GetVoidPointer(0),
                (double*) camDirection->GetVoidPointer(0),
                (double*) camUp->GetVoidPointer(0),
                (double*) camRes->GetVoidPointer(0),
                (double*) camNearFar->GetVoidPointer(0),
                (double*) camHeight->GetVoidPointer(0),
                this->Downsampling,

                // Output
                vertices,
                triangles,
                triangleDistortions
            );
        });
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    vtkSmartPointer<vtkFloatArray> triangleDistortionsScalars = vtkSmartPointer<vtkFloatArray>::New();
    triangleDistortionsScalars->SetNumberOfComponents(1);
    triangleDistortionsScalars->SetName("triangleDistortion");

    {
        for(auto& x: vertices)
            points->InsertNextPoint( get<0>(x), get<1>(x), get<2>(x) );

        mesh->SetPoints(points);
    }

    {
        for(size_t i=0; i<triangles.size(); i++){
            vtkIdType ids[3] = {get<0>(triangles[i]),get<1>(triangles[i]),get<2>(triangles[i])};

            mesh->InsertNextCell(VTK_TRIANGLE, 3, ids);

            triangleDistortionsScalars->InsertTuple1(i, triangleDistortions[i]);
        }

        mesh->GetCellData()->AddArray(triangleDistortionsScalars);
    }

    outputMesh->ShallowCopy(mesh);

    {
        stringstream msg;
        msg << "[ttkDepthImageBasedGeometryApproximation] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
}
