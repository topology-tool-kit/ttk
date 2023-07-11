/// \ingroup base
/// \class ttk::DimensionReduction
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date September 2018
///
/// \brief TTK VTK-filter that apply dimension reduction algorithms on input.
///
/// \param Input Input scalar field (vtkTable)
/// \param Output Output scalar field (vtkTable)
///
/// \brief TTK VTK-filter that takes a matrix (vtkTable) as input and apply a
/// dimension reduction algorithm from scikit-learn.
///
/// \sa ttk::Triangulation
/// \sa ttkDimensionReduction.cpp %for a usage example.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/karhunenLoveDigits64Dimensions//">Karhunen-Love
///   Digits 64-Dimensions example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_householdAnalysis/">Persistent
///   Generators Household Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_periodicPicture/">Persistent
///   Generators Periodic Picture example</a> \n
///

#pragma once

#include <Debug.h>

namespace ttk {

  class DimensionReduction : virtual public Debug {

  public:
    DimensionReduction();

    /** Scikit-Learn Dimension Reduction algorithms */
    enum class METHOD {
      /** Spectral Embedding */
      SE = 0,
      /** Locally Linear Embedding */
      LLE = 1,
      /** Multi-Dimensional Scaling */
      MDS = 2,
      /** t-distributed Stochastic Neighbor Embedding */
      T_SNE = 3,
      /** IsoMap Embedding */
      ISOMAP = 4,
      /** Principal Component Analysis */
      PCA = 5,
    };

    inline void setSEParameters(const std::string &Affinity,
                                const float Gamma,
                                const std::string &EigenSolver,
                                const bool InputIsADistanceMatrix) {
      if(InputIsADistanceMatrix) {
        se_Affinity = "precomputed";
      } else {
        se_Affinity = Affinity;
      }
      se_Gamma = Gamma;
      se_EigenSolver = EigenSolver;
    }

    inline void setLLEParameters(const float Regularization,
                                 const std::string &EigenSolver,
                                 const float Tolerance,
                                 const int MaxIteration,
                                 const std::string &Method_s,
                                 const float HessianTolerance,
                                 const float ModifiedTolerance,
                                 const std::string &NeighborsAlgorithm) {
      lle_Regularization = Regularization;
      lle_EigenSolver = EigenSolver;
      lle_Tolerance = Tolerance;
      lle_MaxIteration = MaxIteration;
      lle_Method = Method_s;
      lle_HessianTolerance = HessianTolerance;
      lle_ModifiedTolerance = ModifiedTolerance;
      lle_NeighborsAlgorithm = NeighborsAlgorithm;
    }

    inline void setMDSParameters(const bool Metric,
                                 const int Init,
                                 const int MaxIteration,
                                 const int Verbose,
                                 const float Epsilon,
                                 const bool Dissimilarity) {
      mds_Metric = Metric;
      mds_Init = Init;
      mds_MaxIteration = MaxIteration;
      mds_Verbose = Verbose;
      mds_Epsilon = Epsilon;
      if(Dissimilarity) {
        mds_Dissimilarity = "precomputed";
      } else {
        mds_Dissimilarity = "euclidean";
      }
    }

    inline void setTSNEParameters(const float Perplexity,
                                  const float Exaggeration,
                                  const float LearningRate,
                                  const int MaxIteration,
                                  const int MaxIterationProgress,
                                  const float GradientThreshold,
                                  const std::string &Metric,
                                  const std::string &Init,
                                  const int Verbose,
                                  const std::string &Method_s,
                                  const float Angle) {
      tsne_Perplexity = Perplexity;
      tsne_Exaggeration = Exaggeration;
      tsne_LearningRate = LearningRate;
      tsne_MaxIteration = MaxIteration;
      tsne_MaxIterationProgress = MaxIterationProgress;
      tsne_GradientThreshold = GradientThreshold;
      tsne_Metric = Metric;
      tsne_Init = Init;
      tsne_Verbose = Verbose;
      tsne_Method = Method_s;
      tsne_Angle = Angle;
    }

    inline void setISOParameters(const std::string &EigenSolver,
                                 const float Tolerance,
                                 const int MaxIteration,
                                 const std::string &PathMethod,
                                 const std::string &NeighborsAlgorithm) {
      iso_EigenSolver = EigenSolver;
      iso_Tolerance = Tolerance;
      iso_MaxIteration = MaxIteration;
      iso_PathMethod = PathMethod;
      iso_NeighborsAlgorithm = NeighborsAlgorithm;
    }

    inline void setPCAParameters(const bool Copy,
                                 const bool Whiten,
                                 const std::string &SVDSolver,
                                 const float Tolerance,
                                 const std::string &MaxIteration) {
      pca_Copy = Copy;
      pca_Whiten = Whiten;
      pca_SVDSolver = SVDSolver;
      pca_Tolerance = Tolerance;
      pca_MaxIteration = MaxIteration;
    }

    inline void setInputModulePath(const std::string &modulePath) {
      ModulePath = modulePath;
    }

    inline void setInputModuleName(const std::string &moduleName) {
      ModuleName = moduleName;
    }

    inline void setInputFunctionName(const std::string &functionName) {
      FunctionName = functionName;
    }

    inline void setInputMethod(METHOD method) {
      this->Method = method;
    }

    inline void setInputNumberOfComponents(const int numberOfComponents) {
      this->NumberOfComponents = numberOfComponents;
    }

    inline void setInputNumberOfNeighbors(const int numberOfNeighbors) {
      this->NumberOfNeighbors = numberOfNeighbors;
    }

    inline void setInputIsDeterministic(const int isDeterm) {
      this->IsDeterministic = isDeterm;
    }

    inline void setIsInputDistanceMatrix(const bool data) {
      if(data) {
        this->se_Affinity = "precomputed";
        this->mds_Dissimilarity = "precomputed";
        this->tsne_Metric = "precomputed";
        this->iso_Metric = "precomputed";
      } else {
        this->se_Affinity = "nearest_neighbors";
        this->mds_Dissimilarity = "euclidean";
        this->tsne_Metric = "euclidean";
        this->iso_Metric = "euclidean";
      }
    }

    bool isPythonFound() const;

    int execute(std::vector<std::vector<double>> &outputEmbedding,
                const std::vector<double> &inputMatrix,
                const int nRows,
                const int nColumns) const;

  protected:
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
    std::string mds_Dissimilarity{"euclidean"};

    // tsne
    float tsne_Perplexity{30};
    float tsne_Exaggeration{12};
    float tsne_LearningRate{200};
    int tsne_MaxIteration{1000};
    int tsne_MaxIterationProgress{300};
    float tsne_GradientThreshold{1e-7};
    std::string tsne_Metric{"euclidean"};
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
    std::string iso_Metric{"euclidean"};

    // pca
    bool pca_Copy{true};
    bool pca_Whiten{false};
    std::string pca_SVDSolver{"auto"};
    float pca_Tolerance{0};
    std::string pca_MaxIteration{"auto"};

    // testing
    std::string ModulePath{"default"};
    std::string ModuleName{"dimensionReduction"};
    std::string FunctionName{"doIt"};

    METHOD Method{METHOD::MDS};
    int NumberOfComponents{2};
    int NumberOfNeighbors{5};
    int IsDeterministic{true};
    char majorVersion_{'0'};
  };
} // namespace ttk
