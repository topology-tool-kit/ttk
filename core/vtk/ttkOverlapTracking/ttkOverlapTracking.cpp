#include <ttkOverlapTracking.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDelimitedTextReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkOverlapTracking)

struct Comparator {
    float* coords;

    inline bool operator() (const size_t& i, const size_t& j) {
        size_t ic = i*3;
        size_t jc = j*3;
        return coords[ic]==coords[jc]
            ? coords[ic+1]==coords[jc+1]
                ? coords[ic+2]<coords[jc+2]
                : coords[ic+1]<coords[jc+1]
            : coords[ic]<coords[jc];
    }
};

int ttkOverlapTracking::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkOverlapTracking] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    // Get Input Object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* inMB = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Get Output Graph
    // vtkInformation* outInfo = outputVector->GetInformationObject(0);
    // vtkUnstructuredGrid* outTable = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* outMB = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    outMB->ShallowCopy(inMB);

    vtkPointSet* block = vtkPointSet::SafeDownCast( inMB->GetBlock(0) );

    if(block==nullptr)
        cout<< "shit"<<endl;

    auto pd = block->GetPointData();

    vtkSmartPointer<vtkDoubleArray> indices = vtkSmartPointer<vtkDoubleArray>::New();
    indices->SetNumberOfComponents(1);
    indices->SetNumberOfValues( block->GetNumberOfPoints() );
    indices->SetName( "Indices" );
    pd->AddArray( indices );

    vtkPoints* points = block->GetPoints();


    vector<size_t> temp( block->GetNumberOfPoints() );
    for(size_t i=0; i<block->GetNumberOfPoints(); i++)
        temp[i] = i;

    Comparator c = Comparator();
    float* pointCoords = (float*) points->GetVoidPointer(0);
    c.coords = pointCoords;
    sort(temp.begin(), temp.end(), c);

    for(size_t i=0; i<block->GetNumberOfPoints(); i++)
        indices->SetValue(temp[i], (double)i);

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkOverlapTracking] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
