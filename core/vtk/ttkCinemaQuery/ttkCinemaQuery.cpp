#include                <ttkCinemaQuery.h>
#include                <vtkStringArray.h>
#include                <vtkVariantArray.h>
#include                <vtkSmartPointer.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaQuery)

int ttkCinemaQuery::doIt(vtkTable* input, vtkTable* output){

    auto* column =  vtkStringArray::SafeDownCast( input->GetColumn(0) );

    string serverAddress = column->GetValue(0);

    vector<vector<string>> rMatrix = cinemaQuery_.execute<int>( serverAddress, QueryString );

    unsigned int rows = rMatrix.size();
    unsigned int cols = rMatrix[0].size();

    if(rows<1)
        return 0;

    for ( unsigned int i = 0; i < cols; i++ ) {
        vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();

        column->SetNumberOfTuples(rows-1);
        column->SetName( rMatrix[0][i].data() );

        for (unsigned int j = 1; j < rows; j++) {
            column->SetValue(j-1, vtkVariant(rMatrix[j][i]) );
        }
        output->AddColumn(column);
    }

    return 0;
}

int ttkCinemaQuery::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkTable* inputTable = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkTable* outputTable = vtkTable::GetData(outputVector);

  doIt(
      inputTable,
      outputTable
  );

  return 1;
}