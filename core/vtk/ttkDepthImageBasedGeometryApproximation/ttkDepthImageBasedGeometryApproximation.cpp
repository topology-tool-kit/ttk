#include                  <ttkDepthImageBasedGeometryApproximation.h>

using namespace std;
using namespace ttk;

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

vtkStandardNewMacro(ttkDepthImageBasedGeometryApproximation)

int ttkDepthImageBasedGeometryApproximation::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkDepthImageBasedGeometryApproximation] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;
    Timer t;

    // Prepare input and output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto inputMBD = vtkMultiBlockDataSet::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outputMBD = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // Set Wrapper
    depthImageBasedGeometryApproximation_.setWrapper(this);

    // Process each depth image individually
    size_t n = inputMBD->GetNumberOfBlocks();
    for(size_t i=0; i<n; i++){
        // Get vtkImageData
        auto inputImage = vtkImageData::SafeDownCast( inputMBD->GetBlock(i) );

        // Get input paramters
        auto depthValues = inputImage->GetPointData()->GetAbstractArray("DepthValues");
        auto camHeight = inputImage->GetFieldData()->GetAbstractArray("CamHeight");
        auto camPosition = inputImage->GetFieldData()->GetAbstractArray("CamPosition");
        auto camDirection = inputImage->GetFieldData()->GetAbstractArray("CamDirection");
        auto camUp = inputImage->GetFieldData()->GetAbstractArray("CamUp");
        auto camNearFar = inputImage->GetFieldData()->GetAbstractArray("CamNearFar");
        auto camRes = inputImage->GetFieldData()->GetAbstractArray("CamRes");

        // Check if all parameters are present
        if(depthValues==nullptr || camHeight==nullptr || camPosition==nullptr || camDirection==nullptr || camUp==nullptr || camNearFar==nullptr || camRes==nullptr){
            stringstream msg;
            msg << "[ttkDepthImageBasedGeometryApproximation] ERROR: Input depth image does not have one or more of the required fields (see Cinema Spec D - Data Product Specification)" << endl;
            dMsg(cout, msg.str(), memoryMsg);
            return 0;
        }

        // Vectors to hold raw approximated geometry
        vector<tuple<double,double,double>> vertices;
        vector<tuple<int,int,int>> triangles;
        vector<double> triangleDistortions;

        // Approximate geometry
        cout<<depthValues->GetDataType()<<endl;
        switch(depthValues->GetDataType()){
            vtkTemplateMacro({
                depthImageBasedGeometryApproximation_.execute<VTK_TT>(
                    // Input
                    (VTK_TT*) depthValues->GetVoidPointer(0),
                    (double*) camPosition->GetVoidPointer(0),
                    (double*) camDirection->GetVoidPointer(0),
                    (double*) camUp->GetVoidPointer(0),
                    (double*) camRes->GetVoidPointer(0),
                    (double*) camNearFar->GetVoidPointer(0),
                    (double*) camHeight->GetVoidPointer(0),
                    this->SubSampling,

                    // Output
                    vertices,
                    triangles,
                    triangleDistortions
                );
            });
        }

        // Represent approximated geometry via VTK
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkSmartPointer<vtkDoubleArray> triangleDistortionsScalars = vtkSmartPointer<vtkDoubleArray>::New();
        triangleDistortionsScalars->SetNumberOfComponents(1);
        triangleDistortionsScalars->SetName("TriangleDistortion");

        // Create points
        {
            for(auto& x: vertices)
                points->InsertNextPoint( get<0>(x), get<1>(x), get<2>(x) );
            mesh->SetPoints(points);
        }

        // Create cells
        {
            for(size_t i=0; i<triangles.size(); i++){
                vtkIdType ids[3] = {get<0>(triangles[i]),get<1>(triangles[i]),get<2>(triangles[i])};
                mesh->InsertNextCell(VTK_TRIANGLE, 3, ids);
                triangleDistortionsScalars->InsertTuple1(i, triangleDistortions[i]);
            }

            mesh->GetCellData()->AddArray(triangleDistortionsScalars);
        }

        outputMBD->SetBlock(i,mesh);
    }

    // Print status
    {
        stringstream msg;
        msg << "[ttkDepthImageBasedGeometryApproximation] ------------------------------------" << endl
            << "[ttkDepthImageBasedGeometryApproximation] " << n << " Images processed" << endl
            << "[ttkDepthImageBasedGeometryApproximation]   time: " << t.getElapsedTime() << " s" << endl
            << "[ttkDepthImageBasedGeometryApproximation] memory: " << m.getElapsedUsage() << " MB" << endl
            << "[ttkDepthImageBasedGeometryApproximation] ------------------------------------" << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
