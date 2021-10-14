#include <ttkContinuousScatterPlot.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCharArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkContinuousScatterPlot);

ttkContinuousScatterPlot::ttkContinuousScatterPlot() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

ttkContinuousScatterPlot::~ttkContinuousScatterPlot() {
}

int ttkContinuousScatterPlot::FillInputPortInformation(int port,
                                                       vtkInformation *info) {

  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  return 1;
}

int ttkContinuousScatterPlot::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {

  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 1;
}

template <typename dataType1, class triangulationType>
int ttkContinuousScatterPlot::dispatch(const dataType1 *scalars1,
                                       vtkDataArray *inputScalars2,
                                       const triangulationType *triangulation) {
  int status = 0;
  switch(inputScalars2->GetDataType()) {
    vtkTemplateMacro(
      (status = this->execute<dataType1, VTK_TT, triangulationType>(
         scalars1, (VTK_TT *)ttkUtils::GetVoidPointer(inputScalars2),
         triangulation)));
  };
  return status;
}

int ttkContinuousScatterPlot::RequestData(vtkInformation *ttkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0], 0);
  vtkDataSet *output = vtkDataSet::GetData(outputVector, 0);
  if(!input || !output)
    return 0;

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  auto inputScalars1 = this->GetInputArrayToProcess(0, input);
  auto inputScalars2 = this->GetInputArrayToProcess(1, input);

#ifndef TTK_ENABLE_KAMIKAZE
  // wrong scalar fields
  if(!inputScalars1 || !inputScalars2) {
    this->printErr("wrong scalar fields.");
    return -2;
  }
#endif

  SimplexId numberOfPixels
    = ScatterplotResolution[0] * ScatterplotResolution[1];
#ifndef TTK_ENABLE_KAMIKAZE
  // no pixels
  if(!numberOfPixels) {
    this->printErr("no pixels.");
    return -3;
  }
#endif

  std::vector<std::vector<double>> density(ScatterplotResolution[0]);
  std::vector<std::vector<char>> validPointMask(ScatterplotResolution[0]);
  for(SimplexId k = 0; k < ScatterplotResolution[0]; ++k) {
    density[k].resize(ScatterplotResolution[1], 0.0);
    validPointMask[k].resize(ScatterplotResolution[1], 0);
  }

  SimplexId numberOfPoints = input->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  // no points
  if(numberOfPoints < 1) {
    this->printErr("no points.");
    return -4;
  }
#endif

  double scalarMin[2];
  double scalarMax[2];
  for(SimplexId k = 0; k < numberOfPoints; ++k) {
    double d1 = inputScalars1->GetTuple1(k);
    double d2 = inputScalars2->GetTuple1(k);

    if(!k or scalarMin[0] > d1)
      scalarMin[0] = d1;
    if(!k or scalarMin[1] > d2)
      scalarMin[1] = d2;
    if(!k or scalarMax[0] < d1)
      scalarMax[0] = d1;
    if(!k or scalarMax[1] < d2)
      scalarMax[1] = d2;
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // scalar fields stats problem
  if(scalarMin[0] == scalarMax[0] or scalarMin[1] == scalarMax[1]) {
    this->printErr("scalar fields stats problem.");
    return -5;
  }
#endif

  // calling the executing package
  this->setVertexNumber(numberOfPoints);
  this->setDummyValue(WithDummyValue, DummyValue);
  this->setResolutions(ScatterplotResolution[0], ScatterplotResolution[1]);
  this->setScalarMin(scalarMin);
  this->setScalarMax(scalarMax);
  this->setOutputDensity(&density);
  this->setOutputMask(&validPointMask);

  int status = 0;
  ttkVtkTemplateMacro(inputScalars1->GetDataType(), triangulation->getType(),
                      (status = this->dispatch<VTK_TT, TTK_TT>(
                         (VTK_TT *)ttkUtils::GetVoidPointer(inputScalars1),
                         inputScalars2, (TTK_TT *)triangulation->getData())));

  // something wrong in baseCode
  if(status != 0) {
    std::stringstream msg;
    msg << "ContinuousScatterPlot.execute() error " << status;
    this->printErr(msg.str());
    return -6;
  }

  vtkNew<vtkCharArray> maskScalars;
  maskScalars->SetNumberOfComponents(1);
  maskScalars->SetNumberOfTuples(numberOfPixels);
  maskScalars->SetName("ValidPointMask");

  vtkNew<vtkDoubleArray> densityScalars;
  densityScalars->SetNumberOfComponents(1);
  densityScalars->SetNumberOfTuples(numberOfPixels);
  densityScalars->SetName("Density");

  vtkNew<vtkDoubleArray> scalars1;
  scalars1->SetNumberOfComponents(1);
  scalars1->SetNumberOfTuples(numberOfPixels);
  scalars1->SetName(inputScalars1->GetName());

  vtkNew<vtkDoubleArray> scalars2;
  scalars2->SetNumberOfComponents(1);
  scalars2->SetNumberOfTuples(numberOfPixels);
  scalars2->SetName(inputScalars2->GetName());

  double delta[2];
  delta[0] = (scalarMax[0] - scalarMin[0]) / (ScatterplotResolution[0] - 1);
  delta[1] = (scalarMax[1] - scalarMin[1]) / (ScatterplotResolution[1] - 1);

  double imageMin[2]{0.0, 0.0};
  double imageMax[2]{1.0, 1.0};
  if(ProjectImageSupport) {
    imageMin[0] = scalarMin[0];
    imageMin[1] = scalarMin[1];
    imageMax[0] = scalarMax[0];
    imageMax[1] = scalarMax[1];
  }
  double imageDelta[2];
  imageDelta[0] = (imageMax[0] - imageMin[0]) / (ScatterplotResolution[0] - 1);
  imageDelta[1] = (imageMax[1] - imageMin[1]) / (ScatterplotResolution[1] - 1);

  vtkNew<vtkUnstructuredGrid> vtu;
  vtkNew<vtkPoints> pts;
  pts->SetNumberOfPoints(numberOfPixels);

  SimplexId id{};
  vtkIdType ids[3];
  for(SimplexId i = 0; i < ScatterplotResolution[0]; i++) {
    for(SimplexId j = 0; j < ScatterplotResolution[1]; j++) {
      // positions:
      double x = imageMin[0] + i * imageDelta[0];
      double y = imageMin[1] + j * imageDelta[1];
      pts->SetPoint(id, x, y, 0);

      // scalars:
      // valid mask
      maskScalars->SetTuple1(id, validPointMask[i][j]);
      // density
      densityScalars->SetTuple1(id, density[i][j]);
      // original scalar fields
      double d1 = scalarMin[0] + i * delta[0];
      double d2 = scalarMin[1] + j * delta[1];
      scalars1->SetTuple1(id, d1);
      scalars2->SetTuple1(id, d2);
      if(i < ScatterplotResolution[0] - 1
         and j < ScatterplotResolution[1] - 1) {
        ids[0] = id;
        ids[1] = id + 1;
        ids[2] = id + ScatterplotResolution[1];
        vtu->InsertNextCell(VTK_TRIANGLE, 3, ids);

        ids[0] = id + 1;
        ids[1] = id + ScatterplotResolution[1];
        ids[2] = id + ScatterplotResolution[1] + 1;
        vtu->InsertNextCell(VTK_TRIANGLE, 3, ids);
      }

      ++id;
    }
  }
  vtu->SetPoints(pts);
  vtu->GetPointData()->AddArray(maskScalars);
  vtu->GetPointData()->AddArray(densityScalars);
  vtu->GetPointData()->AddArray(scalars1);
  vtu->GetPointData()->AddArray(scalars2);
  output->ShallowCopy(vtu);

  return 1;
}
