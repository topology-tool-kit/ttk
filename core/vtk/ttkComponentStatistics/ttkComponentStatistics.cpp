#include                  <ttkComponentStatistics.h>
#include                  <vtkPointData.h>
#include                  <vtkFieldData.h>
#include                  <vtkMultiBlockDataSet.h>
#include                  <vtkDataSet.h>
#include                  <vtkDoubleArray.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkComponentStatistics)

int ttkComponentStatistics::RequestData (
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    Timer t;
    Memory m;

    // Print Status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkComponentStatistics] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* inputMBD = vtkMultiBlockDataSet::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    auto block = vtkDataSet::SafeDownCast( inputMBD->GetBlock(0) );
    auto pointData = block->GetPointData();

    auto vertexNumbers = pointData->GetArray("VertexNumber");
    auto regionIDs = pointData->GetArray("RegionId");


    double* regionLimits = regionIDs->GetRange();

    cout<<regionLimits[0]<<" "<<regionLimits[1]<<endl;
    cout<<vertexNumbers->GetSize()<<endl;
    cout<<regionIDs->GetSize()<<endl;

    size_t n = vertexNumbers->GetSize();
    double* vertexNumbersData = (double*) vertexNumbers->GetVoidPointer(0);
    vtkIdType* regionIDsData = (vtkIdType*) regionIDs->GetVoidPointer(0);


    vtkSmartPointer<vtkDoubleArray> v = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> r = vtkSmartPointer<vtkDoubleArray>::New();
    v->SetName("VertexNumber");
    r->SetName("RegionID");
    v->SetNumberOfValues(regionLimits[1]+1);
    r->SetNumberOfValues(regionLimits[1]+1);

    for(size_t i=0; i<r->GetSize(); i++){
        r->SetValue(i, (double)i);
    }

    for(size_t i=0; i<n; i++){
        v->SetValue(regionIDsData[i], vertexNumbersData[i]);
    }

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkTable* outTable = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkTable> newTable = vtkSmartPointer<vtkTable>::New();
    newTable->AddColumn( r );
    newTable->AddColumn( v );
    outTable->ShallowCopy(newTable);

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkComponentStatistics]   Time: " << t.getElapsedTime() << " s." << endl;
        msg << "[ttkComponentStatistics] Memory: " << m.getElapsedUsage() << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
