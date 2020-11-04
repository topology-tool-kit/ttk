#include <ttkCinemaDarkroomColorMapping.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

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

int ttkCinemaDarkroomColorMapping::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomColorMapping::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template<typename DT>
int mapScalarsToColor(
  unsigned char* color,
  const std::vector<double>& colorMap,
  const double* nanColor,
  const DT* array,
  const size_t nPixels
){
  const size_t nKeys = colorMap.size()/4;

  for(size_t i=0; i<nPixels; i++){

    const double value = (double)array[i];
    if(isnan(value)){
      size_t idx = i*3;
      color[idx+0] = 255.0*nanColor[0];
      color[idx+1] = 255.0*nanColor[1];
      color[idx+2] = 255.0*nanColor[2];
      continue;
    }

    size_t ki = 0;
    for(size_t k=1; k<nKeys; k++){
      if(value<=colorMap[k*4]){
        ki = k-1;
        break;
      }
    }

    double lambda = (value-colorMap[ki*4]) / (colorMap[(ki+1)*4]-colorMap[ki*4]);
    double lambdaInv = 1-lambda;

    size_t idx = i*3;
    size_t idx2 = ki*4;
    color[idx+0] = 255.0 * (lambdaInv*colorMap[idx2+1] + lambda*colorMap[idx2+5]);
    color[idx+1] = 255.0 * (lambdaInv*colorMap[idx2+2] + lambda*colorMap[idx2+6]);
    color[idx+2] = 255.0 * (lambdaInv*colorMap[idx2+3] + lambda*colorMap[idx2+7]);
  }
  return 1;
};


int ttkCinemaDarkroomColorMapping::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  ttk::Timer timer;
  this->printMsg("Applying Color Map",0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);

  auto input = vtkImageData::GetData(inputVector[0]);
  auto output = vtkImageData::GetData(outputVector);
  output->ShallowCopy(input);

  auto scalarArray = this->GetInputArrayToProcess(0, output);
  if(!scalarArray || this->GetInputArrayAssociation(0,output)!=0 || scalarArray->GetNumberOfComponents()!=1){
    this->printErr("Unable to retrieve point scalar array.");
    return 0;
  }

  size_t nPixels = scalarArray->GetNumberOfTuples();

  auto colorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colorArray->SetName( "Diffuse" );
  colorArray->SetNumberOfComponents(3);
  colorArray->SetNumberOfTuples( nPixels );
  output->GetPointData()->AddArray( colorArray );
  auto colorArrayData = static_cast<unsigned char*>(ttkUtils::GetVoidPointer(colorArray));


  std::vector<double> manualColorMap;
  const std::vector<double>* colorMap{nullptr};

  if(this->ColorMap==-1){ // Solid
    manualColorMap.resize(8);
    manualColorMap[0] = 0;
    manualColorMap[1] = this->SolidColor[0];
    manualColorMap[2] = this->SolidColor[1];
    manualColorMap[3] = this->SolidColor[2];
    manualColorMap[4] = 1;
    manualColorMap[5] = this->SolidColor[0];
    manualColorMap[6] = this->SolidColor[1];
    manualColorMap[7] = this->SolidColor[2];
    colorMap = &manualColorMap;
  } else if(this->ColorMap==-2) {
    int status = ttkUtils::stringListToDoubleVector(this->ManualColorMap, manualColorMap);
    if(!status || manualColorMap.size()<8 || manualColorMap.size()%4!=0){
      this->printErr("Invalid manual color map input.");
      return 0;
    }
    colorMap = &manualColorMap;
  } else {
    if(this->ColorMap<0 || this->ColorMap>=(int)this->ColorMaps.size()){
      this->printErr("Invalid color map index: "+std::to_string(this->ColorMap));
      return 0;
    }
    colorMap = &this->ColorMaps[this->ColorMap];
  }

  switch(scalarArray->GetDataType()){
    vtkTemplateMacro(
      mapScalarsToColor<VTK_TT>(
        colorArrayData,
        *colorMap,
        this->NANColor,
        static_cast<VTK_TT*>(ttkUtils::GetVoidPointer(scalarArray)),
        nPixels
      )
    );
  }

  this->printMsg("Applying Color Map",1,timer.getElapsedTime(),this->threadNumber_);

  return 1;
}