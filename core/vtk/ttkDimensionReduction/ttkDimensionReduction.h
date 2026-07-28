/// \class ttkDimensionReduction
/// \ingroup vtk
/// \author GuillaumeFavelier <guillaume.favelier@gmail.com>
/// \date September 2018.
///
/// \brief TTK VTK-filter that wraps the ttk::DimensionReduction
/// processing package.
///
/// VTK wrapping code for the ttk::DimensionReduction package.
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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/karhunenLoveDigits64Dimensions/">Karhunen-Love
///   Digits 64-Dimensions example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramPGA/">Persistence
///   Diagram Principal Geodesic Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_householdAnalysis/">Persistent
///   Generators Household Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_periodicPicture/">Persistent
///   Generators Periodic Picture example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topoAEppTeaser/">Topological
///   Autoencoders++ Teaser example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topoMapTeaser/">TopoMap
///   Teaser example</a> \n
///

/// \b Related \b publication: \n
/// "Topomap: A 0-dimensional homology preserving projection of high-dimensional
/// data"\n Harish Doraiswamy, Julien Tierny, Paulo J. S. Silva, Luis Gustavo
/// Nonato, and Claudio Silva\n Proc. of IEEE VIS 2020.\n IEEE Transactions on
/// Visualization and Computer Graphics 27(2): 561-571, 2020. \n
///
/// "Topological Autoencoders" \n
/// Michael Moor, Max Horn, Bastian Rieck, Karsten Borgwardt, \n
/// Proceedings of the 37th International Conference on Machine Learning,
/// 2020. \n
///
/// "Optimizing persistent homology-based functions" \n
/// Mathieu Carriere, Frederic Chazal, Marc Glisse, Yuichi Ike,
/// Hariprasad Kannan, Yuhei Umeda, \n
/// Proceedings of the 38th International Conference on Machine Learning,
/// 2021. \n
///
/// "Topological Autoencoders++: Fast and Accurate Cycle-Aware Dimensionality
/// Reduction" \n
/// MattÃ©o ClÃ©mot, Julie Digne, Julien Tierny, \n
/// IEEE Transactions on Visualization and Computer Graphics.
/// Accepted, to be presented at IEEE VIS 2026.

#pragma once

// VTK Module
#include <ttkDimensionReductionModule.h>

// TTK includes
#include <DimensionReduction.h>
#include <TopoMap.h>
#include <TopologicalDimensionReduction.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class TTKDIMENSIONREDUCTION_EXPORT ttkDimensionReduction
  : public ttkAlgorithm,
    protected ttk::DimensionReduction {

public:
  static ttkDimensionReduction *New();
  vtkTypeMacro(ttkDimensionReduction, ttkAlgorithm);

  void SetScalarFields(const std::string &s) {
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

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

  void SetInitializationFields(const std::string &s) {
    InitializationFields.push_back(s);
    Modified();
  }

  void ClearInitializationFields() {
    InitializationFields.clear();
    Modified();
  }

  vtkSetMacro(SelectInitializationFieldsWithRegexp, bool);
  vtkGetMacro(SelectInitializationFieldsWithRegexp, bool);

  vtkSetMacro(InitializationRegexpString, const std::string &);
  vtkGetMacro(InitializationRegexpString, std::string);

  vtkSetMacro(NumberOfComponents, int);
  vtkGetMacro(NumberOfComponents, int);

  vtkSetMacro(NumberOfNeighbors, int);
  vtkGetMacro(NumberOfNeighbors, int);

  vtkSetMacro(IsDeterministic, int);
  vtkGetMacro(IsDeterministic, int);

  ttkSetEnumMacro(Method, METHOD);
  vtkGetEnumMacro(Method, METHOD);

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

  // SE && MDS
  void SetInputIsADistanceMatrix(const bool b) {
    this->InputIsADistanceMatrix = b;
    this->setIsInputDistanceMatrix(b);
    Modified();
  }
  vtkGetMacro(InputIsADistanceMatrix, bool);

  // SE
  vtkSetMacro(se_Affinity, const std::string &);
  vtkGetMacro(se_Affinity, std::string);

  vtkSetMacro(se_Gamma, float);
  vtkGetMacro(se_Gamma, float);

  vtkSetMacro(se_EigenSolver, const std::string &);
  vtkGetMacro(se_EigenSolver, std::string);

  // LLE
  vtkSetMacro(lle_Regularization, float);
  vtkGetMacro(lle_Regularization, float);

  vtkSetMacro(lle_EigenSolver, const std::string &);
  vtkGetMacro(lle_EigenSolver, std::string);

  vtkSetMacro(lle_Tolerance, float);
  vtkGetMacro(lle_Tolerance, float);

  vtkSetMacro(lle_MaxIteration, int);
  vtkGetMacro(lle_MaxIteration, int);

  vtkSetMacro(lle_Method, const std::string &);
  vtkGetMacro(lle_Method, std::string);

  vtkSetMacro(lle_HessianTolerance, float);
  vtkGetMacro(lle_HessianTolerance, float);

  vtkSetMacro(lle_ModifiedTolerance, float);
  vtkGetMacro(lle_ModifiedTolerance, float);

  vtkSetMacro(lle_NeighborsAlgorithm, const std::string &);
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

  vtkSetMacro(tsne_Metric, const std::string &);
  vtkGetMacro(tsne_Metric, std::string);

  vtkSetMacro(tsne_Init, const std::string &);
  vtkGetMacro(tsne_Init, std::string);

  vtkSetMacro(tsne_Verbose, int);
  vtkGetMacro(tsne_Verbose, int);

  vtkSetMacro(tsne_Method, const std::string &);
  vtkGetMacro(tsne_Method, std::string);

  vtkSetMacro(tsne_Angle, float);
  vtkGetMacro(tsne_Angle, float);

  // Iso
  vtkSetMacro(iso_EigenSolver, const std::string &);
  vtkGetMacro(iso_EigenSolver, std::string);

  vtkSetMacro(iso_Tolerance, float);
  vtkGetMacro(iso_Tolerance, float);

  vtkSetMacro(iso_MaxIteration, int);
  vtkGetMacro(iso_MaxIteration, int);

  vtkSetMacro(iso_PathMethod, const std::string &);
  vtkGetMacro(iso_PathMethod, std::string);

  vtkSetMacro(iso_NeighborsAlgorithm, const std::string &);
  vtkGetMacro(iso_NeighborsAlgorithm, std::string);

  vtkSetMacro(iso_Metric, const std::string &);
  vtkGetMacro(iso_Metric, std::string);

  // PCA
  vtkSetMacro(pca_Copy, bool);
  vtkGetMacro(pca_Copy, bool);

  vtkSetMacro(pca_Whiten, bool);
  vtkGetMacro(pca_Whiten, bool);

  vtkSetMacro(pca_SVDSolver, const std::string &);
  vtkGetMacro(pca_SVDSolver, std::string);

  vtkSetMacro(pca_Tolerance, float);
  vtkGetMacro(pca_Tolerance, float);

  vtkSetMacro(pca_MaxIteration, const std::string &);
  vtkGetMacro(pca_MaxIteration, std::string);

  // TopoMap
  vtkSetMacro(topomap_AngularSampleNb, unsigned long int);
  vtkGetMacro(topomap_AngularSampleNb, unsigned long int);

  vtkSetMacro(topomap_CheckMST, bool);
  vtkGetMacro(topomap_CheckMST, bool);

  ttkSetEnumMacro(topomap_Strategy, ttk::TopoMap::STRATEGY);
  vtkGetEnumMacro(topomap_Strategy, ttk::TopoMap::STRATEGY);

  // AutoEncoder
  vtkSetMacro(ae_CUDA, bool);
  vtkGetMacro(ae_CUDA, bool);

  vtkSetMacro(ae_Deterministic, bool);
  vtkGetMacro(ae_Deterministic, bool);

  vtkSetMacro(ae_Seed, int);
  vtkGetMacro(ae_Seed, int);

  vtkSetMacro(ae_Epochs, int);
  vtkGetMacro(ae_Epochs, int);

  vtkSetMacro(ae_LearningRate, double);
  vtkGetMacro(ae_LearningRate, double);

  ttkSetEnumMacro(ae_Method, ttk::TopologicalDimensionReduction::REGUL);
  vtkGetEnumMacro(ae_Method, ttk::TopologicalDimensionReduction::REGUL);

  ttkSetEnumMacro(ae_Optimizer, ttk::TopologicalDimensionReduction::OPTIMIZER);
  vtkGetEnumMacro(ae_Optimizer, ttk::TopologicalDimensionReduction::OPTIMIZER);

  ttkSetEnumMacro(ae_Model, ttk::TopologicalDimensionReduction::MODEL);
  vtkGetEnumMacro(ae_Model, ttk::TopologicalDimensionReduction::MODEL);

  vtkSetMacro(ae_Architecture, const std::string &);
  vtkGetMacro(ae_Architecture, std::string);

  vtkSetMacro(ae_Activation, const std::string &);
  vtkGetMacro(ae_Activation, std::string);

  vtkSetMacro(ae_BatchSize, int);
  vtkGetMacro(ae_BatchSize, int);

  vtkSetMacro(ae_BatchNormalization, bool);
  vtkGetMacro(ae_BatchNormalization, bool);

  vtkSetMacro(ae_RegCoefficient, double);
  vtkGetMacro(ae_RegCoefficient, double);

  vtkSetMacro(IsInputImages, bool);
  vtkGetMacro(IsInputImages, bool);

  vtkSetMacro(ae_PreOptimize, bool);
  vtkGetMacro(ae_PreOptimize, bool);

  vtkSetMacro(ae_PreOptimizeEpochs, int);
  vtkGetMacro(ae_PreOptimizeEpochs, int);

  // testing
  vtkSetMacro(ModulePath, const std::string &);
  vtkGetMacro(ModulePath, std::string);

  vtkSetMacro(ModuleName, const std::string &);
  vtkGetMacro(ModuleName, std::string);

  vtkSetMacro(FunctionName, const std::string &);
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
  std::vector<std::string> ScalarFields{};

  bool SelectInitializationFieldsWithRegexp{false};
  std::string InitializationRegexpString{".*"};
  std::vector<std::string> InitializationFields{};

  bool KeepAllDataArrays{true};

  // mds && se
  bool InputIsADistanceMatrix{false};

  std::vector<std::vector<double>> outputData_{};
};
