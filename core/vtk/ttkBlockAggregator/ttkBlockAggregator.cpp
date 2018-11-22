#include <ttkBlockAggregator.h>

#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiBlockDataSet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkBlockAggregator)

int ttkBlockAggregator::AggregateBlock(vtkDataObject* dataObject, vtkMultiBlockDataSet* mb, size_t index, bool useShallowCopy){
    auto temp0 = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    temp0->SetBlock(0, dataObject);

    auto temp1 = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    if(useShallowCopy)
        temp1->ShallowCopy( temp0 );
    else
        temp1->DeepCopy( temp0 );

    mb->SetBlock(index, temp1->GetBlock(0));

    return 1;
}

int ttkBlockAggregator::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg << "================================================================================" << endl
            << "[ttkBlockAggregator] RequestData" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Get Input and Output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

    // Get iteration information
    double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );

    // First timestep
    if(this->GetForceReset() || i==0 || this->AggregatedMultiBlockDataSet==nullptr)
        this->AggregatedMultiBlockDataSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    size_t nBlocks = this->AggregatedMultiBlockDataSet->GetNumberOfBlocks();
    for(size_t i=0; i<5; i++){
        vtkInformationVector* inVector = inputVector[i];

        if(inVector->GetNumberOfInformationObjects()!=1) continue;


        vtkInformation* inInfo = inVector->GetInformationObject(0);
        auto input = inInfo->Get( vtkDataObject::DATA_OBJECT() );
        if(input){
            size_t index = nBlocks + i;

            // Add block
            this->AggregateBlock(input, this->AggregatedMultiBlockDataSet, index, true);

            // Print status
            stringstream msg;
            msg << "[ttkBlockAggregator] Added block at index " << index << endl;
            dMsg(cout, msg.str(), infoMsg);
        }
    }

    // Get Output
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outMB = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );
    outMB->ShallowCopy( this->AggregatedMultiBlockDataSet );

    return 1;
}