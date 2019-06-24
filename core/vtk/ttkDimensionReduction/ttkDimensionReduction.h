/// \ingroup vtk
/// \class ttkDimensionReduction
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the dimensionReduction processing package.
///
/// VTK wrapping code for the @DimensionReduction package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::DimensionReduction
#pragma once

#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkTableAlgorithm.h>

#include <DimensionReduction.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDimensionReduction
#else
class ttkDimensionReduction
#endif
  : public vtkTableAlgorithm,
    public ttk::Wrapper {
public:
  enum Method {
    SpectralEmbedding = 0,
    LocallyLinearEmbedding,
    MDS,
    TSNE,
    Isomap,
    PCA
  };

  static ttkDimensionReduction *New();
  vtkTypeMacro(ttkDimensionReduction, vtkTableAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  void SetScalarFields(std::string s) {
    ScalarFields.push_back(s);
    Modified();
  }

  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  // default
  vtkSetMacro(NumberOfComponents, int);
  vtkGetMacro(NumberOfComponents, int);

  vtkSetMacro(NumberOfNeighbors, int);
  vtkGetMacro(NumberOfNeighbors, int);

  vtkSetMacro(IsDeterministic, int);
  vtkGetMacro(IsDeterministic, int);

  vtkSetMacro(Method, int);
  vtkGetMacro(Method, int);

  vtkSetMacro(KeepAllDataArrays, int);
  vtkGetMacro(KeepAllDataArrays, int);

  // SE
  vtkSetMacro(se_Affinity, std::string);
  vtkGetMacro(se_Affinity, std::string);

  vtkSetMacro(se_Gamma, float);
  vtkGetMacro(se_Gamma, float);

  vtkSetMacro(se_EigenSolver, std::string);
  vtkGetMacro(se_EigenSolver, std::string);

  // LLE
  vtkSetMacro(lle_Regularization, float);
  vtkGetMacro(lle_Regularization, float);

  vtkSetMacro(lle_EigenSolver, std::string);
  vtkGetMacro(lle_EigenSolver, std::string);

  vtkSetMacro(lle_Tolerance, float);
  vtkGetMacro(lle_Tolerance, float);

  vtkSetMacro(lle_MaxIteration, int);
  vtkGetMacro(lle_MaxIteration, int);

  vtkSetMacro(lle_Method, std::string);
  vtkGetMacro(lle_Method, std::string);

  vtkSetMacro(lle_HessianTolerance, float);
  vtkGetMacro(lle_HessianTolerance, float);

  vtkSetMacro(lle_ModifiedTolerance, float);
  vtkGetMacro(lle_ModifiedTolerance, float);

  vtkSetMacro(lle_NeighborsAlgorithm, std::string);
  vtkGetMacro(lle_NeighborsAlgorithm, std::string);

  // MDS
  vtkSetMacro(mds_Metric, bool);
  vtkGetMacro(mds_Metric, bool);

  vtkSetMacro(mds_Init, int);
  vtkGetMacro(mds_Init, int);

  vtkSetMacro(mds_MaxIteration, int);
  vtkGetMacro(mds_MaxIteration, int);

  vtkSetMacro(mds_Verbose, int);
  vtkGetMacro(mds_Verbose, int);

  vtkSetMacro(mds_Epsilon, float);
  vtkGetMacro(mds_Epsilon, float);

  vtkSetMacro(mds_Dissimilarity, std::string);
  vtkGetMacro(mds_Dissimilarity, std::string);

  // TSNE
  vtkSetMacro(tsne_Perplexity, float);
  vtkGetMacro(tsne_Perplexity, float);

  vtkSetMacro(tsne_Exaggeration, float);
  vtkGetMacro(tsne_Exaggeration, float);

  vtkSetMacro(tsne_LearningRate, float);
  vtkGetMacro(tsne_LearningRate, float);

  vtkSetMacro(tsne_MaxIteration, int);
  vtkGetMacro(tsne_MaxIteration, int);

  vtkSetMacro(tsne_MaxIterationProgress, int);
  vtkGetMacro(tsne_MaxIterationProgress, int);

  vtkSetMacro(tsne_GradientThreshold, float);
  vtkGetMacro(tsne_GradientThreshold, float);

  vtkSetMacro(tsne_Metric, std::string);
  vtkGetMacro(tsne_Metric, std::string);

  vtkSetMacro(tsne_Init, std::string);
  vtkGetMacro(tsne_Init, std::string);

  vtkSetMacro(tsne_Verbose, int);
  vtkGetMacro(tsne_Verbose, int);

  vtkSetMacro(tsne_Method, std::string);
  vtkGetMacro(tsne_Method, std::string);

  vtkSetMacro(tsne_Angle, float);
  vtkGetMacro(tsne_Angle, float);

  // Iso
  vtkSetMacro(iso_EigenSolver, std::string);
  vtkGetMacro(iso_EigenSolver, std::string);

  vtkSetMacro(iso_Tolerance, float);
  vtkGetMacro(iso_Tolerance, float);

  vtkSetMacro(iso_MaxIteration, int);
  vtkGetMacro(iso_MaxIteration, int);

  vtkSetMacro(iso_PathMethod, std::string);
  vtkGetMacro(iso_PathMethod, std::string);

  vtkSetMacro(iso_NeighborsAlgorithm, std::string);
  vtkGetMacro(iso_NeighborsAlgorithm, std::string);

  // PCA
  vtkSetMacro(pca_Copy, bool);
  vtkGetMacro(pca_Copy, bool);

  vtkSetMacro(pca_Whiten, bool);
  vtkGetMacro(pca_Whiten, bool);

  vtkSetMacro(pca_SVDSolver, std::string);
  vtkGetMacro(pca_SVDSolver, std::string);

  vtkSetMacro(pca_Tolerance, float);
  vtkGetMacro(pca_Tolerance, float);

  vtkSetMacro(pca_MaxIteration, std::string);
  vtkGetMacro(pca_MaxIteration, std::string);

  // testing
  vtkSetMacro(ModulePath, std::string);
  vtkGetMacro(ModulePath, std::string);

  vtkSetMacro(ModuleName, std::string);
  vtkGetMacro(ModuleName, std::string);

  vtkSetMacro(FunctionName, std::string);
  vtkGetMacro(FunctionName, std::string);

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
        break;
    }

    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
        break;
    }

    return 1;
  }

protected:
  ttkDimensionReduction() {
    NumberOfComponents = 2;
    NumberOfNeighbors = 5;
    Method = 2;

    se_Affinity = "nearest_neighbors";
    se_Gamma = 1;
    se_EigenSolver = "auto";

    lle_Regularization = 1e-3;
    lle_EigenSolver = "auto";
    lle_Tolerance = 1e-3;
    lle_MaxIteration = 300;
    lle_Method = "standard";
    lle_HessianTolerance = 1e-3;
    lle_ModifiedTolerance = 1e-3;
    lle_NeighborsAlgorithm = "auto";

    mds_Metric = true;
    mds_Init = 4;
    mds_MaxIteration = 300;
    mds_Verbose = 0;
    mds_Epsilon = 0;
    mds_Dissimilarity = "euclidean";

    tsne_Perplexity = 30;
    tsne_Exaggeration = 12;
    tsne_LearningRate = 200;
    tsne_MaxIteration = 1000;
    tsne_MaxIterationProgress = 300;
    tsne_GradientThreshold = 1e-7;
    tsne_Metric = "euclidean";
    tsne_Init = "random";
    tsne_Verbose = 0;
    tsne_Method = "barnes_hut";
    tsne_Angle = 0.5;

    iso_EigenSolver = "auto";
    iso_Tolerance = 1e-3;
    iso_MaxIteration = 300;
    iso_PathMethod = "auto";
    iso_NeighborsAlgorithm = "auto";

    pca_Copy = true;
    pca_Whiten = false;
    pca_SVDSolver = "auto";
    pca_Tolerance = 0;
    pca_MaxIteration = "auto";

    UseAllCores = true;
  }

  ~ttkDimensionReduction() {
  }

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int doIt(vtkTable *input, vtkTable *output);
  bool needsToAbort() override;
  int updateProgress(const float &progress) override;

  // default
  int NumberOfComponents;
  int NumberOfNeighbors;
  int Method;
  int IsDeterministic;
  bool KeepAllDataArrays;

  // se
  std::string se_Affinity;
  float se_Gamma;
  std::string se_EigenSolver;

  // lle
  float lle_Regularization;
  std::string lle_EigenSolver;
  float lle_Tolerance;
  int lle_MaxIteration;
  std::string lle_Method;
  float lle_HessianTolerance;
  float lle_ModifiedTolerance;
  std::string lle_NeighborsAlgorithm;

  // mds
  bool mds_Metric;
  int mds_Init;
  int mds_MaxIteration;
  int mds_Verbose;
  float mds_Epsilon;
  std::string mds_Dissimilarity;

  // tsne
  float tsne_Perplexity;
  float tsne_Exaggeration;
  float tsne_LearningRate;
  int tsne_MaxIteration;
  int tsne_MaxIterationProgress;
  float tsne_GradientThreshold;
  std::string tsne_Metric;
  std::string tsne_Init;
  int tsne_Verbose;
  std::string tsne_Method;
  float tsne_Angle;

  // iso
  std::string iso_EigenSolver;
  float iso_Tolerance;
  int iso_MaxIteration;
  std::string iso_PathMethod;
  std::string iso_NeighborsAlgorithm;

  // pca
  bool pca_Copy;
  bool pca_Whiten;
  std::string pca_SVDSolver;
  float pca_Tolerance;
  std::string pca_MaxIteration;

  // testing
  std::string ModulePath;
  std::string ModuleName;
  std::string FunctionName;
  bool UseAllCores;
  ttk::ThreadId ThreadNumber;
  ttk::DimensionReduction dimensionReduction_;

  std::vector<std::string> ScalarFields;
  std::vector<std::vector<double>> outputData_{};
};
