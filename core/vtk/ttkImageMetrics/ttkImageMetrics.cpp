#include  <ttkImageMetrics.h>

#include  <vtkMultiBlockDataSet.h>
#include  <vtkDataSet.h>
#include  <vtkPointData.h>
#include  <vtkTable.h>
#include  <vtkDoubleArray.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkImageMetrics)

int ttkImageMetrics::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkImageMetrics] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    // Get ImagesA
    vtkInformation* inInfoA = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* imagesA = vtkMultiBlockDataSet::SafeDownCast(inInfoA->Get(vtkDataObject::DATA_OBJECT()));

    // Get ImagesB
    vtkInformation* inInfoB = inputVector[1]->GetInformationObject(0);
    vtkMultiBlockDataSet* imagesB = vtkMultiBlockDataSet::SafeDownCast(inInfoB->Get(vtkDataObject::DATA_OBJECT()));

    // Get Output Table
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkTable* outTable = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Compute Metrics
    {
        size_t n = imagesA->GetNumberOfBlocks();

        // ADD
        {
            vtkSmartPointer<vtkDoubleArray> adds = vtkSmartPointer<vtkDoubleArray>::New();
            adds->SetNumberOfComponents(1);
            adds->SetNumberOfValues(n);
            adds->SetName("ADD");

            for(size_t i=0; i<n; i++){
                auto a = vtkDataSet::SafeDownCast( imagesA->GetBlock(i) );
                auto b = vtkDataSet::SafeDownCast( imagesB->GetBlock(i) );
                auto pixelsA = a->GetPointData()->GetAbstractArray("DepthValues");
                auto pixelsB = b->GetPointData()->GetAbstractArray("DepthValues");

                double add;

                switch(pixelsA->GetDataType()){
                    case VTK_FLOAT: // Float
                        imageMetrics.ADD<float>(
                            // Input
                            (float*) pixelsA->GetVoidPointer(0),
                            (float*) pixelsB->GetVoidPointer(0),
                            pixelsA->GetSize(),

                            // Output
                            add
                        );
                        break;
                    case VTK_DOUBLE: // Double
                        imageMetrics.ADD<double>(
                            // Input
                            (double*) pixelsA->GetVoidPointer(0),
                            (double*) pixelsB->GetVoidPointer(0),
                            pixelsA->GetSize(),

                            // Output
                            add
                        );
                        break;
                    default:
                        // Print Error
                        {
                            stringstream msg;
                            msg << "[ttkImageMetrics] ERROR: Unsorported pixel datatype for ADD-Metric"
                                << m.getElapsedUsage()
                                << " MB." << endl;
                            dMsg(cout, msg.str(), memoryMsg);
                        }
                }

                adds->SetValue(i, add);
            }

            outTable->AddColumn( adds );
        }
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkImageMetrics] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}