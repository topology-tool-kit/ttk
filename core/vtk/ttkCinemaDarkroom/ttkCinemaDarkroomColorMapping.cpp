#include "BaseClass.h"
#include <ttkCinemaDarkroomColorMapping.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <math.h>

vtkStandardNewMacro(ttkCinemaDarkroomColorMapping);

ttkCinemaDarkroomColorMapping::ttkCinemaDarkroomColorMapping() {
  this->setDebugMsgPrefix("CinemaDarkroomColorMapping");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomColorMapping::~ttkCinemaDarkroomColorMapping() {
}

template <typename DT>
int mapScalarsToColor(unsigned char *color,
                      const std::vector<double> &colorMap,
                      const double *nanColor,
                      const DT *array,
                      const double range[2],
                      const size_t nPixels,
                      const int threadNumber) {
  const size_t nKeys = colorMap.size() / 4;
  const double valueDelta = range[1] - range[0];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif
  for(size_t i = 0; i < nPixels; i++) {

    const double value = (double)array[i];
    if(isnan(value)) {
      size_t idx = i * 3;
      color[idx + 0] = 255.0 * nanColor[0];
      color[idx + 1] = 255.0 * nanColor[1];
      color[idx + 2] = 255.0 * nanColor[2];
      continue;
    }

    const double normalizedValue
      = std::max(0.0, std::min(1.0, (value - range[0]) / valueDelta));

    size_t ki = 0;
    for(size_t k = 1; k < nKeys; k++) {
      if(normalizedValue <= colorMap[k * 4]) {
        ki = k - 1;
        break;
      }
    }

    double lambda = (normalizedValue - colorMap[ki * 4])
                    / (colorMap[(ki + 1) * 4] - colorMap[ki * 4]);
    double lambdaInv = 1 - lambda;

    size_t idx = i * 3;
    size_t idx2 = ki * 4;
    color[idx + 0]
      = 255.0 * (lambdaInv * colorMap[idx2 + 1] + lambda * colorMap[idx2 + 5]);
    color[idx + 1]
      = 255.0 * (lambdaInv * colorMap[idx2 + 2] + lambda * colorMap[idx2 + 6]);
    color[idx + 2]
      = 255.0 * (lambdaInv * colorMap[idx2 + 3] + lambda * colorMap[idx2 + 7]);
  }

  TTK_FORCE_USE(threadNumber);
  return 1;
}

int ttkCinemaDarkroomColorMapping::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  ttk::Timer timer;
  this->printMsg("Applying Color Map", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  auto input = vtkImageData::GetData(inputVector[0]);
  auto output = vtkImageData::GetData(outputVector);
  output->ShallowCopy(input);

  auto scalarArray = this->GetInputArrayToProcess(0, output);
  if(!scalarArray || this->GetInputArrayAssociation(0, output) != 0
     || scalarArray->GetNumberOfComponents() != 1) {
    this->printErr("Unable to retrieve point scalar array.");
    return 0;
  }

  size_t nPixels = scalarArray->GetNumberOfTuples();

  double range[2];
  scalarArray->GetRange(range);

  auto colorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colorArray->SetName("Diffuse");
  colorArray->SetNumberOfComponents(3);
  colorArray->SetNumberOfTuples(nPixels);
  output->GetPointData()->AddArray(colorArray);
  auto colorArrayData
    = static_cast<unsigned char *>(ttkUtils::GetVoidPointer(colorArray));

  std::vector<double> manualColorMap;
  const std::vector<double> *colorMap{nullptr};

  if(this->ColorMap == -1) { // Solid
    manualColorMap.resize(8);
    manualColorMap[0] = 0;
    manualColorMap[1] = this->SingleColor[0];
    manualColorMap[2] = this->SingleColor[1];
    manualColorMap[3] = this->SingleColor[2];
    manualColorMap[4] = 1;
    manualColorMap[5] = this->SingleColor[0];
    manualColorMap[6] = this->SingleColor[1];
    manualColorMap[7] = this->SingleColor[2];
    colorMap = &manualColorMap;
  } else if(this->ColorMap == -2) {
    int status = ttkUtils::stringListToDoubleVector(
      this->ManualColorMap, manualColorMap);
    if(!status || manualColorMap.size() < 8 || manualColorMap.size() % 4 != 0) {
      this->printErr("Invalid manual color map input.");
      return 0;
    }
    colorMap = &manualColorMap;
  } else {
    if(this->ColorMap < 0 || this->ColorMap >= (int)this->ColorMaps.size()) {
      this->printErr("Invalid color map index: "
                     + std::to_string(this->ColorMap));
      return 0;
    }
    colorMap = &this->ColorMaps[this->ColorMap];
  }

  switch(scalarArray->GetDataType()) {
    vtkTemplateMacro(mapScalarsToColor<VTK_TT>(
      colorArrayData, *colorMap, this->NANColor,
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)), range,
      nPixels, this->threadNumber_));
  }

  this->printMsg(
    "Applying Color Map", 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}
