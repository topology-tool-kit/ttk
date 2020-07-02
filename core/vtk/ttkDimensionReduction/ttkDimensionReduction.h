/// \class ttkDimensionReduction
/// \ingroup vtk
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

// VTK Module
#include <ttkDimensionReductionModule.h>

// TTK includes
#include <DimensionReduction.h>
#include <ttkAlgorithm.h>

class TTKDIMENSIONREDUCTION_EXPORT ttkDimensionReduction
  : public ttkAlgorithm,
    protected ttk::DimensionReduction {

public:
  static ttkDimensionReduction *New();
  vtkTypeMacro(ttkDimensionReduction, ttkAlgorithm);

  void SetScalarFields(std::string s) {
    ScalarFields.push_back(s);
    Modified();
  }

  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  // default
  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(RegexpString, std::string);
  vtkGetMacro(RegexpString, std::string);

  vtkSetMacro(NumberOfComponents, int);
  vtkGetMacro(NumberOfComponents, int);

  vtkSetMacro(NumberOfNeighbors, int);
  vtkGetMacro(NumberOfNeighbors, int);

  vtkSetMacro(IsDeterministic, int);
  vtkGetMacro(IsDeterministic, int);

  vtkSetMacro(Method, int);
  vtkGetMacro(Method, int);

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

  // SE && MDS
  vtkSetMacro(InputIsADistanceMatrix, bool);
  vtkGetMacro(InputIsADistanceMatrix, bool);

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

protected:
  ttkDimensionReduction();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // default
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  int NumberOfComponents{2};
  int NumberOfNeighbors{5};
  int Method{2}; // MDS
  int IsDeterministic{true};
  bool KeepAllDataArrays{true};

  // mds && se
  bool InputIsADistanceMatrix{false};

  // se
  std::string se_Affinity{"nearest_neighbors"};
  float se_Gamma{1};
  std::string se_EigenSolver{"auto"};

  // lle
  float lle_Regularization{1e-3};
  std::string lle_EigenSolver{"auto"};
  float lle_Tolerance{1e-3};
  int lle_MaxIteration{300};
  std::string lle_Method{"standard"};
  float lle_HessianTolerance{1e-3};
  float lle_ModifiedTolerance{1e-3};
  std::string lle_NeighborsAlgorithm{"auto"};

  // mds
  bool mds_Metric{true};
  int mds_Init{4};
  int mds_MaxIteration{300};
  int mds_Verbose{0};
  float mds_Epsilon{0};

  // tsne
  float tsne_Perplexity{30};
  float tsne_Exaggeration{12};
  float tsne_LearningRate{200};
  int tsne_MaxIteration{1000};
  int tsne_MaxIterationProgress{300};
  float tsne_GradientThreshold{1e-7};
  std::string tsne_Metric{"euclidian"};
  std::string tsne_Init{"random"};
  int tsne_Verbose{0};
  std::string tsne_Method{"barnes_hut"};
  float tsne_Angle{0.5};

  // iso
  std::string iso_EigenSolver{"auto"};
  float iso_Tolerance{1e-3};
  int iso_MaxIteration{300};
  std::string iso_PathMethod{"auto"};
  std::string iso_NeighborsAlgorithm{"auto"};

  // pca
  bool pca_Copy{true};
  bool pca_Whiten{false};
  std::string pca_SVDSolver{"auto"};
  float pca_Tolerance{0};
  std::string pca_MaxIteration{"auto"};

  // testing
  std::string ModulePath{};
  std::string ModuleName{};
  std::string FunctionName{};

  std::vector<std::string> ScalarFields{};
  std::vector<std::vector<double>> outputData_{};
};
