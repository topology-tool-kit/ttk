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

#pragma once

#include <Debug.h>

namespace ttk {

  class DimensionReduction : virtual public Debug {

  public:
    DimensionReduction();

    inline int setSEParameters(std::string &Affinity,
                               float Gamma,
                               std::string &EigenSolver,
                               bool InputIsADistanceMatrix) {
      if(InputIsADistanceMatrix) {
        se_Affinity = "precomputed";
      } else {
        se_Affinity = Affinity;
      }
      se_Gamma = Gamma;
      se_EigenSolver = EigenSolver;
      return 0;
    }

    inline int setLLEParameters(float Regularization,
                                std::string &EigenSolver,
                                float Tolerance,
                                int MaxIteration,
                                std::string &Method,
                                float HessianTolerance,
                                float ModifiedTolerance,
                                std::string &NeighborsAlgorithm) {
      lle_Regularization = Regularization;
      lle_EigenSolver = EigenSolver;
      lle_Tolerance = Tolerance;
      lle_MaxIteration = MaxIteration;
      lle_Method = Method;
      lle_HessianTolerance = HessianTolerance;
      lle_ModifiedTolerance = ModifiedTolerance;
      lle_NeighborsAlgorithm = NeighborsAlgorithm;
      return 0;
    }

    inline int setMDSParameters(bool Metric,
                                int Init,
                                int MaxIteration,
                                int Verbose,
                                float Epsilon,
                                bool Dissimilarity) {
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
      return 0;
    }

    inline int setTSNEParameters(float Perplexity,
                                 float Exaggeration,
                                 float LearningRate,
                                 int MaxIteration,
                                 int MaxIterationProgress,
                                 float GradientThreshold,
                                 std::string &Metric,
                                 std::string &Init,
                                 int Verbose,
                                 std::string &Method,
                                 float Angle) {
      tsne_Perplexity = Perplexity;
      tsne_Exaggeration = Exaggeration;
      tsne_LearningRate = LearningRate;
      tsne_MaxIteration = MaxIteration;
      tsne_MaxIterationProgress = MaxIterationProgress;
      tsne_GradientThreshold = GradientThreshold;
      tsne_Metric = Metric;
      tsne_Init = Init;
      tsne_Verbose = Verbose;
      tsne_Method = Method;
      tsne_Angle = Angle;
      return 0;
    }

    inline int setISOParameters(std::string &EigenSolver,
                                float Tolerance,
                                int MaxIteration,
                                std::string &PathMethod,
                                std::string &NeighborsAlgorithm) {
      iso_EigenSolver = EigenSolver;
      iso_Tolerance = Tolerance;
      iso_MaxIteration = MaxIteration;
      iso_PathMethod = PathMethod;
      iso_NeighborsAlgorithm = NeighborsAlgorithm;
      return 0;
    }

    inline int setPCAParameters(bool Copy,
                                bool Whiten,
                                std::string &SVDSolver,
                                float Tolerance,
                                std::string &MaxIteration) {
      pca_Copy = Copy;
      pca_Whiten = Whiten;
      pca_SVDSolver = SVDSolver;
      pca_Tolerance = Tolerance;
      pca_MaxIteration = MaxIteration;
      return 0;
    }

    inline int setInputModulePath(const std::string &modulePath) {
      ModulePath = modulePath;
      return 0;
    }

    inline int setInputModuleName(const std::string &moduleName) {
      ModuleName = moduleName;
      return 0;
    }

    inline int setInputFunctionName(const std::string &functionName) {
      FunctionName = functionName;
      return 0;
    }

    inline int setInputMatrixDimensions(SimplexId numberOfRows,
                                        SimplexId numberOfColumns) {
      numberOfRows_ = numberOfRows;
      numberOfColumns_ = numberOfColumns;
      return 0;
    }

    inline int setInputMatrix(void *data) {
      matrix_ = data;
      return 0;
    }

    inline int setInputMethod(int method) {
      method_ = method;
      return 0;
    }

    inline int setInputNumberOfComponents(int numberOfComponents) {
      numberOfComponents_ = numberOfComponents;
      return 0;
    }

    inline int setInputNumberOfNeighbors(int numberOfNeighbors) {
      numberOfNeighbors_ = numberOfNeighbors;
      return 0;
    }

    inline int setInputIsDeterministic(int randomState) {
      randomState_ = randomState;
      return 0;
    }

    inline int setOutputComponents(std::vector<std::vector<double>> *data) {
      embedding_ = data;
      return 0;
    }

    bool isPythonFound() const;

    int execute() const;

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
    std::string ModulePath{};
    std::string ModuleName{};
    std::string FunctionName{};

    SimplexId numberOfRows_{0};
    SimplexId numberOfColumns_{0};
    int method_{};
    int numberOfComponents_{0};
    int numberOfNeighbors_{0};
    int randomState_{0};
    void *matrix_{};
    std::vector<std::vector<double>> *embedding_{};
    char majorVersion_{'0'};
  };
} // namespace ttk
