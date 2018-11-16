#include <ttkDepthImageBasedGeometryApproximation.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <vtkDataArray.h>
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
        msg<<"================================================================================"<<endl;
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
        auto inPointData = inputImage->GetPointData();
        auto depthValues = inputImage->GetPointData()->GetAbstractArray( this->GetDepthScalarField().data() );

        auto inFieldData = inputImage->GetFieldData();
        auto camHeight = inFieldData->GetAbstractArray("CamHeight");
        auto camPosition = inFieldData->GetAbstractArray("CamPosition");
        auto camDirection = inFieldData->GetAbstractArray("CamDirection");
        auto camUp = inFieldData->GetAbstractArray("CamUp");
        auto camNearFar = inFieldData->GetAbstractArray("CamNearFar");
        auto camRes = inFieldData->GetAbstractArray("CamRes");

        // Check if all parameters are present
        if(depthValues==nullptr || camHeight==nullptr || camPosition==nullptr || camDirection==nullptr || camUp==nullptr || camNearFar==nullptr || camRes==nullptr){
            stringstream msg;
            msg << "[ttkDepthImageBasedGeometryApproximation] ERROR: Input depth image does not have one or more of the required fields (see Cinema Spec D - Data Product Specification)" << endl;
            dMsg(cout, msg.str(), memoryMsg);
            return 0;
        }

        // Vectors to hold raw approximated geometry
        vector<size_t> indicies;
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
                    this->GetSubsampling(),

                    // Output
                    indicies,
                    vertices,
                    triangles,
                    triangleDistortions
                );
            });
        }

        // Represent approximated geometry via VTK
        auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

        // Create points
        {
            size_t n = vertices.size();

            auto points = vtkSmartPointer<vtkPoints>::New();
            points->SetNumberOfPoints( n );

            auto pointCoords = (float*) points->GetVoidPointer(0);
            size_t i=0;
            for(auto& x: vertices){
                pointCoords[i++] = get<0>(x);
                pointCoords[i++] = get<1>(x);
                pointCoords[i++] = get<2>(x);
            }

            mesh->SetPoints( points );
        }

        // Copy Point Data
        {
            size_t n = inPointData->GetNumberOfArrays();

            auto outPointData = mesh->GetPointData();
            size_t m = indicies.size();

            for(size_t i=0; i<n; i++){
                auto inArray = inPointData->GetArray(i);

                auto outArray = vtkDataArray::CreateDataArray( inArray->GetDataType() );
                outArray->SetName( inArray->GetName() );
                outArray->SetNumberOfTuples( m );

                for(size_t j=0; j<m; j++){
                    outArray->SetTuple(j, inArray->GetTuple(indicies[j]));
                }

                outPointData->AddArray( outArray );
            }
        }

        // Create cells
        {
            size_t n = triangles.size();

            auto cells = vtkSmartPointer<vtkIdTypeArray>::New();
            cells->SetNumberOfValues( 4*n );
            auto cellIds = (vtkIdType*) cells->GetVoidPointer(0);

            auto triangleDistortionsScalars = vtkSmartPointer<vtkDoubleArray>::New();
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

            auto cellArray = vtkSmartPointer<vtkCellArray>::New();
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
        msg << "[ttkDepthImageBasedGeometryApproximation] --------------------------------------" << endl
            << "[ttkDepthImageBasedGeometryApproximation] " << n << " Images processed" << endl
            << "[ttkDepthImageBasedGeometryApproximation]   time: " << t.getElapsedTime() << " s" << endl
            << "[ttkDepthImageBasedGeometryApproximation] memory: " << m.getElapsedUsage() << " MB" << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
