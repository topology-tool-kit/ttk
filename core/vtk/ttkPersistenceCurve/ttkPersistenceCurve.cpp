#include <ttkPersistenceCurve.h>

using namespace std;
using namespace ttk;

using namespace ftm;

vtkStandardNewMacro(ttkPersistenceCurve)

  ttkPersistenceCurve::ttkPersistenceCurve()
  : UseAllCores{}, inputScalars_{}, JTPersistenceCurve_{vtkTable::New()},
    MSCPersistenceCurve_{vtkTable::New()}, STPersistenceCurve_{vtkTable::New()},
    CTPersistenceCurve_{vtkTable::New()}, offsets_{}, inputOffsets_{},
    varyingMesh_{} {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(4);

  triangulation_ = NULL;
  ScalarFieldId = 0;
  OffsetFieldId = -1;
  InputOffsetScalarFieldName = ttk::OffsetScalarFieldName;
  ForceInputOffsetScalarField = false;
  ComputeSaddleConnectors = false;
  UseAllCores = true;

  inputTriangulation_ = vtkSmartPointer<ttkTriangulationFilter>::New();
}

ttkPersistenceCurve::~ttkPersistenceCurve() {
  if(JTPersistenceCurve_)
    JTPersistenceCurve_->Delete();
  if(MSCPersistenceCurve_)
    MSCPersistenceCurve_->Delete();
  if(STPersistenceCurve_)
    STPersistenceCurve_->Delete();
  if(CTPersistenceCurve_)
    CTPersistenceCurve_->Delete();
  if(offsets_)
    offsets_->Delete();
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkPersistenceCurve::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPersistenceCurve::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkPersistenceCurve] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

vtkTable *ttkPersistenceCurve::GetOutput() {
  return this->GetOutput(0);
}

vtkTable *ttkPersistenceCurve::GetOutput(int port) {
  return vtkTable::SafeDownCast(this->GetOutputDataObject(port));
}

int ttkPersistenceCurve::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;

    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;

    case 2:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;

    case 3:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;
  }

  return 1;
}

int ttkPersistenceCurve::getScalars(vtkDataSet *input) {
  vtkPointData *pointData = input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    cerr << "[ttkPersistenceCurve] Error : input has no point data." << endl;
    return -1;
  }
#endif

  if(ScalarField.length()) {
    inputScalars_ = pointData->GetArray(ScalarField.data());
  } else {
    inputScalars_ = pointData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_) {
    cerr << "[ttkPersistenceCurve] Error : input scalar field pointer is null."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkPersistenceCurve::getTriangulation(vtkDataSet *input) {

  varyingMesh_ = false;

  triangulation_ = ttkTriangulation::getTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_) {
    cerr << "[ttkPersistenceCurve] Error: input triangulation is NULL." << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  persistenceCurve_.setupTriangulation(triangulation_);

  if(triangulation_->isEmpty()
     or ttkTriangulation::hasChangedConnectivity(triangulation_, input, this)) {
    Modified();
    varyingMesh_ = true;
  }

  return 0;
}

int ttkPersistenceCurve::getOffsets(vtkDataSet *input) {
  if(OffsetFieldId != -1) {
    inputOffsets_ = input->GetPointData()->GetArray(OffsetFieldId);
    if(inputOffsets_) {
      InputOffsetScalarFieldName = inputOffsets_->GetName();
      ForceInputOffsetScalarField = true;
    }
  }

  if(ForceInputOffsetScalarField and InputOffsetScalarFieldName.length()) {
    inputOffsets_
      = input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  } else if(input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    inputOffsets_ = input->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  } else {
    if(varyingMesh_ and offsets_) {
      offsets_->Delete();
      offsets_ = nullptr;
    }

    if(!offsets_) {
      const SimplexId numberOfVertices = input->GetNumberOfPoints();

      offsets_ = ttkSimplexIdTypeArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName(ttk::OffsetScalarFieldName);
      for(SimplexId i = 0; i < numberOfVertices; ++i)
        offsets_->SetTuple1(i, i);
    }

    inputOffsets_ = offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_) {
    cerr << "[ttkPersistenceCurve] Error : wrong input offset scalar field."
         << endl;
    return -1;
  }
#endif

  return 0;
}
int ttkPersistenceCurve::doIt(vtkDataSet *input,
                              vtkTable *outputJTPersistenceCurve,
                              vtkTable *outputMSCPersistenceCurve,
                              vtkTable *outputSTPersistenceCurve,
                              vtkTable *outputCTPersistenceCurve) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkPersistenceCurve] Error: input pointer is NULL." << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkPersistenceCurve] Error: input has no point." << endl;
    return -1;
  }

  if(!outputJTPersistenceCurve or !outputMSCPersistenceCurve
     or !outputSTPersistenceCurve or !outputCTPersistenceCurve) {
    cerr << "[ttkPersistenceCurve] Error: output pointer is NULL." << endl;
    return -1;
  }
#endif

  int ret{};

  ret = getScalars(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkPersistenceCurve] Error : wrong scalars." << endl;
    return -1;
  }
#endif

  ret = getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkPersistenceCurve] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret = getOffsets(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkPersistenceCurve] Error : wrong offsets." << endl;
    return -1;
  }
  if(inputOffsets_->GetDataType() != VTK_INT
     and inputOffsets_->GetDataType() != VTK_ID_TYPE) {
    cerr
      << "[ttkPersistenceCurve] Error : input offset field type not supported."
      << endl;
    return -1;
  }
#endif

  persistenceCurve_.setWrapper(this);
  persistenceCurve_.setInputScalars(inputScalars_->GetVoidPointer(0));
  persistenceCurve_.setInputOffsets(inputOffsets_->GetVoidPointer(0));
  persistenceCurve_.setComputeSaddleConnectors(ComputeSaddleConnectors);
  switch(inputScalars_->GetDataType()) {
#ifndef _MSC_VER
    vtkTemplateMacro(({
      vector<pair<VTK_TT, SimplexId>> JTPlot;
      vector<pair<VTK_TT, SimplexId>> STPlot;
      vector<pair<VTK_TT, SimplexId>> MSCPlot;
      vector<pair<VTK_TT, SimplexId>> CTPlot;

      persistenceCurve_.setOutputJTPlot(&JTPlot);
      persistenceCurve_.setOutputMSCPlot(&MSCPlot);
      persistenceCurve_.setOutputSTPlot(&STPlot);
      persistenceCurve_.setOutputCTPlot(&CTPlot);
      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = persistenceCurve_.execute<VTK_TT, int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = persistenceCurve_.execute<VTK_TT, vtkIdType>();

#ifndef TTK_ENABLE_KAMIKAZE
      if(ret) {
        cerr << "[ttkPersistenceCurve] PersistenceCurve.execute() error code : "
             << ret << endl;
        return -1;
      }
#endif

      ret = getPersistenceCurve<vtkDoubleArray, VTK_TT>(TreeType::Join, JTPlot);
#ifndef TTK_ENABLE_KAMIKAZE
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error :"
             << " build of join tree persistence curve has failed." << endl;
        return -1;
      }
#endif

      ret = getMSCPersistenceCurve<vtkDoubleArray, VTK_TT>(MSCPlot);
#ifndef TTK_ENABLE_KAMIKAZE
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error : "
             << "build of saddle-saddle persistence curve has failed." << endl;
        return -1;
      }
#endif

      ret
        = getPersistenceCurve<vtkDoubleArray, VTK_TT>(TreeType::Split, STPlot);
#ifndef TTK_ENABLE_KAMIKAZE
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error : "
             << "build of split tree persistence curve has failed." << endl;
        return -1;
      }
#endif

      ret = getPersistenceCurve<vtkDoubleArray, VTK_TT>(
        TreeType::Contour, CTPlot);
#ifndef TTK_ENABLE_KAMIKAZE
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error : "
             << "build of contour tree persistence curve has failed." << endl;
        return -1;
      }
#endif
    }));
#else
#ifndef TTK_ENABLE_KAMIKAZE
    vtkTemplateMacro({
      vector<pair<VTK_TT TTK_COMMA SimplexId>> JTPlot;
      vector<pair<VTK_TT TTK_COMMA SimplexId>> STPlot;
      vector<pair<VTK_TT TTK_COMMA SimplexId>> MSCPlot;
      vector<pair<VTK_TT TTK_COMMA SimplexId>> CTPlot;

      persistenceCurve_.setOutputJTPlot(&JTPlot);
      persistenceCurve_.setOutputMSCPlot(&MSCPlot);
      persistenceCurve_.setOutputSTPlot(&STPlot);
      persistenceCurve_.setOutputCTPlot(&CTPlot);
      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = persistenceCurve_.execute<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = persistenceCurve_.execute<VTK_TT TTK_COMMA vtkIdType>();

      if(ret) {
        cerr << "[ttkPersistenceCurve] PersistenceCurve.execute() error code : "
             << ret << endl;
        return -1;
      }

      ret = getPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(
        TreeType::Join, JTPlot);
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error :"
             << " build of join tree persistence curve has failed." << endl;
        return -1;
      }

      ret = getMSCPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(MSCPlot);
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error : "
             << "build of saddle-saddle persistence curve has failed." << endl;
        return -1;
      }

      ret = getPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(
        TreeType::Split, STPlot);
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error : "
             << "build of split tree persistence curve has failed." << endl;
        return -1;
      }

      ret = getPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(
        TreeType::Contour, CTPlot);
      if(ret) {
        cerr << "[ttkPersistenceCurve] Error : "
             << "build of contour tree persistence curve has failed." << endl;
        return -1;
      }
    });
#else
    vtkTemplateMacro({
      vector<pair<VTK_TT TTK_COMMA SimplexId>> JTPlot;
      vector<pair<VTK_TT TTK_COMMA SimplexId>> STPlot;
      vector<pair<VTK_TT TTK_COMMA SimplexId>> MSCPlot;
      vector<pair<VTK_TT TTK_COMMA SimplexId>> CTPlot;

      persistenceCurve_.setOutputJTPlot(&JTPlot);
      persistenceCurve_.setOutputMSCPlot(&MSCPlot);
      persistenceCurve_.setOutputSTPlot(&STPlot);
      persistenceCurve_.setOutputCTPlot(&CTPlot);
      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = persistenceCurve_.execute<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = persistenceCurve_.execute<VTK_TT TTK_COMMA vtkIdType>();

      ret = getPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(
        TreeType::Join, JTPlot);

      ret = getMSCPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(MSCPlot);

      ret = getPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(
        TreeType::Split, STPlot);

      ret = getPersistenceCurve<vtkDoubleArray TTK_COMMA VTK_TT>(
        TreeType::Contour, CTPlot);
    });
#endif
#endif
  }

  outputJTPersistenceCurve->ShallowCopy(JTPersistenceCurve_);
  outputMSCPersistenceCurve->ShallowCopy(MSCPersistenceCurve_);
  outputSTPersistenceCurve->ShallowCopy(STPersistenceCurve_);
  outputCTPersistenceCurve->ShallowCopy(CTPersistenceCurve_);

  return 0;
}

int ttkPersistenceCurve::RequestData(vtkInformation *request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  Memory m;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);

  inputTriangulation_->SetInputData(input);
  inputTriangulation_->Update();

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkTable *outputJTPersistenceCurve
    = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo = outputVector->GetInformationObject(1);
  vtkTable *outputMSCPersistenceCurve
    = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo = outputVector->GetInformationObject(2);
  vtkTable *outputSTPersistenceCurve
    = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo = outputVector->GetInformationObject(3);
  vtkTable *outputCTPersistenceCurve
    = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  doIt(inputTriangulation_->GetOutput(), outputJTPersistenceCurve,
       outputMSCPersistenceCurve, outputSTPersistenceCurve,
       outputCTPersistenceCurve);

  {
    stringstream msg;
    msg << "[ttkPersistenceCurve] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
