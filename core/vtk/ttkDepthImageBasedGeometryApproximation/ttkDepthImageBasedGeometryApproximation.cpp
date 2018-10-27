#include <ttkDepthImageBasedGeometryApproximation.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDepthImageBasedGeometryApproximation)

int ttkDepthImageBasedGeometryApproximation::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    Memory m;
    Timer t;

    // Print status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkDepthImageBasedGeometryApproximation] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

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
        auto depthValues = inputImage->GetPointData()->GetAbstractArray("Depth");
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
                    this->Subsampling,

                    // Output
                    vertices,
                    triangles,
                    triangleDistortions
                );
            });
        }

        // Represent approximated geometry via VTK
        vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

        // Create points
        {
            size_t n = vertices.size();

            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            points->SetNumberOfPoints( n );

            float* pointCoords = (float*) points->GetVoidPointer(0);
            size_t i=0;
            for(auto& x: vertices){
                pointCoords[i++] = get<0>(x);
                pointCoords[i++] = get<1>(x);
                pointCoords[i++] = get<2>(x);
            }

            mesh->SetPoints( points );
        }

        // Create cells
        {
            size_t n = triangles.size();

            vtkSmartPointer<vtkIdTypeArray> cells = vtkSmartPointer<vtkIdTypeArray>::New();
            cells->SetNumberOfValues(n + n * 3);
            vtkIdType* cellIds = (vtkIdType*) cells->GetVoidPointer(0);

            vtkSmartPointer<vtkDoubleArray> triangleDistortionsScalars = vtkSmartPointer<vtkDoubleArray>::New();
            triangleDistortionsScalars->SetNumberOfValues( n );
            triangleDistortionsScalars->SetNumberOfComponents(1);
            triangleDistortionsScalars->SetName("Distortion");
            double* triangleDistortionsScalarsData = (double*) triangleDistortionsScalars->GetVoidPointer(0);

            size_t q=0;
            for(size_t i=0; i<triangles.size(); i++){
                cellIds[q++] = 3;
                cellIds[q++] = (vtkIdType) get<0>(triangles[i]);
                cellIds[q++] = (vtkIdType) get<1>(triangles[i]);
                cellIds[q++] = (vtkIdType) get<2>(triangles[i]);

                triangleDistortionsScalarsData[i] = triangleDistortions[i];
            }

            vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
            cellArray->SetCells(n, cells);
            mesh->SetCells(VTK_TRIANGLE, cellArray);

            mesh->GetCellData()->AddArray(triangleDistortionsScalars);
        }

        outputMBD->SetBlock(i,mesh);
        this->updateProgress( ((float)i)/((float)(n-1)) );
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
