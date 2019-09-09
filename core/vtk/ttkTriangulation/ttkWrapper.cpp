#include <ttkWrapper.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTriangulationFilter) // constructor
  ttkTriangulationFilter::ttkTriangulationFilter() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkTriangulationFilter::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkTriangulationFilter::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkTriangulationFilter] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkTriangulationFilter::RequestDataObject(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSetAlgorithm::RequestDataObject(request, inputVector, outputVector);

  for(int i = 0; i < GetNumberOfOutputPorts(); i++) {

    vtkInformation *outputInformation = outputVector->GetInformationObject(i);

    vtkDataSet *output = vtkDataSet::SafeDownCast(
      outputInformation->Get(vtkDataObject::DATA_OBJECT()));

    if(output) {

      string dataType = output->GetClassName();

      TTK_POLY_DATA_NEW(i, outputInformation, dataType);

      TTK_UNSTRUCTURED_GRID_NEW(i, outputInformation, dataType);

      // NOTE:
      // We do not allocate a ttkImageData here because paraview seems to be
      // confused by them (no PointData or CellData to refer to).
      // this is not a problem since everything is implicit here anyway.
      // (but we do need to allocate a ttkImageData in the input of TTK filters
      // to access the triangulation data-structure).
      if(dataType == "vtkImageData") {
        ttkImageData *data = ttkImageData::SafeDownCast(
          outputInformation->Get(vtkDataObject::DATA_OBJECT()));

        if(!data) {
          data = ttkImageData::New();
          outputInformation->Set(vtkDataObject::DATA_OBJECT(), data);
          data->FastDelete();
          data->CopyInformationFromPipeline(outputInformation);
          GetOutputPortInformation(i)->Set(
            vtkDataObject::DATA_EXTENT_TYPE(), data->GetExtentType());
        }

        if((data) && (!data->getTriangulation())) {
          data->allocate();
        }
      }
    }
  }

  return 1;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkTriangulationFilter::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  output->ShallowCopy(input);

  return 1;
}
