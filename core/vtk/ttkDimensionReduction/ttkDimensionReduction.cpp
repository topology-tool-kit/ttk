#include <ttkDimensionReduction.h>
#include <ttkUtils.h>

#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkDimensionReduction);

ttkDimensionReduction::ttkDimensionReduction() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkDimensionReduction::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDimensionReduction::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDimensionReduction::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  using ttk::SimplexId;

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  vtkTable *output = vtkTable::GetData(outputVector);

  if(SelectFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    ScalarFields.clear();
    const auto n = input->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = input->GetColumnName(i);
      if(std::regex_match(name, std::regex(RegexpString))) {
        ScalarFields.emplace_back(name);
      }
    }
  }

  if(dimensionReduction_.isPythonFound()) {
    const SimplexId numberOfRows = input->GetNumberOfRows();
    const SimplexId numberOfColumns = ScalarFields.size();

#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfRows <= 0 or numberOfColumns <= 0) {
      this->printErr("input matrix has invalid dimensions");
      return -1;
    }
#endif

    std::vector<double> inputData;
    std::vector<vtkAbstractArray *> arrays;
    for(auto s : ScalarFields)
      arrays.push_back(input->GetColumnByName(s.data()));
    for(SimplexId i = 0; i < numberOfRows; ++i) {
      for(auto arr : arrays)
        inputData.push_back(arr->GetVariantValue(i).ToDouble());
    }

    outputData_.clear();

    dimensionReduction_.setInputModulePath(ModulePath);
    dimensionReduction_.setInputModuleName(ModuleName);
    dimensionReduction_.setInputFunctionName(FunctionName);
    dimensionReduction_.setInputMatrixDimensions(numberOfRows, numberOfColumns);
    dimensionReduction_.setInputMatrix(inputData.data());
    dimensionReduction_.setInputMethod(Method);
    dimensionReduction_.setInputNumberOfComponents(NumberOfComponents);
    dimensionReduction_.setInputNumberOfNeighbors(NumberOfNeighbors);
    dimensionReduction_.setInputIsDeterministic(IsDeterministic);
    dimensionReduction_.setSEParameters(
      se_Affinity, se_Gamma, se_EigenSolver, InputIsADistanceMatrix);
    dimensionReduction_.setLLEParameters(
      lle_Regularization, lle_EigenSolver, lle_Tolerance, lle_MaxIteration,
      lle_Method, lle_HessianTolerance, lle_ModifiedTolerance,
      lle_NeighborsAlgorithm);
    dimensionReduction_.setMDSParameters(mds_Metric, mds_Init, mds_MaxIteration,
                                         mds_Verbose, mds_Epsilon,
                                         InputIsADistanceMatrix);
    dimensionReduction_.setTSNEParameters(
      tsne_Perplexity, tsne_Exaggeration, tsne_LearningRate, tsne_MaxIteration,
      tsne_MaxIterationProgress, tsne_GradientThreshold, tsne_Metric, tsne_Init,
      tsne_Verbose, tsne_Method, tsne_Angle);
    dimensionReduction_.setISOParameters(iso_EigenSolver, iso_Tolerance,
                                         iso_MaxIteration, iso_PathMethod,
                                         iso_NeighborsAlgorithm);
    dimensionReduction_.setPCAParameters(
      pca_Copy, pca_Whiten, pca_SVDSolver, pca_Tolerance, pca_MaxIteration);
    dimensionReduction_.setOutputComponents(&outputData_);
    const int errorCode = dimensionReduction_.execute();

    if(!errorCode) {
      if(KeepAllDataArrays)
        output->ShallowCopy(input);

      for(int i = 0; i < NumberOfComponents; ++i) {
        std::string s = "Component_" + std::to_string(i);
        vtkNew<vtkDoubleArray> arr{};
        ttkUtils::SetVoidArray(arr, outputData_[i].data(), numberOfRows, 1);
        arr->SetName(s.data());
        output->AddColumn(arr);
      }
    }
  } else {
    output->ShallowCopy(input);
    this->printWrn("Python/Numpy not found, features are disabled.");
    vtkWarningMacro("[ttkDimensionReduction] Warning: Python/Numpy not found, "
                    "features are disabled.");
  }

  return 1;
}
