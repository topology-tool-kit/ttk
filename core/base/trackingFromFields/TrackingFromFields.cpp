#include <TrackingFromFields.h>

ttk::TrackingFromFields::TrackingFromFields()
{
  inputData_ = nullptr;
  numberOfInputs_ = 0;
}

ttk::TrackingFromFields::~TrackingFromFields()
{
  if (inputData_)
    free(inputData_);
}

int ttk::TrackingFromFields::performDiagramComputation()
{
//  for (int i = 0; i < (int) fieldNumber; ++i)
//  {
//    ttkPersistenceDiagram *ttkPD = ttkPersistenceDiagram::New();
//
//    ttkPD->AddInputData(input);
//    ttkPD->SetUseInputOffsetScalarField(false);
//    ttkPD->SetComputeSaddleConnectors(false);
//    ttkPD->SetShowInsideDomain(true);
//    ttkPD->SetThreadNumber(1);
//    // ttkPD->SetScalarFieldId(i);
//    ttkPD->SetScalarField(inputScalarFields[i]->GetName());
//    ttkPD->Update();
//
//    vtkUnstructuredGrid *pd =
//        vtkUnstructuredGrid::SafeDownCast(ttkPD->GetOutput(0));
//
//    persistenceDiagrams[i] = pd;
//  }

  ttk::PersistenceDiagram persistenceDiagram_;
  return 0;
}
