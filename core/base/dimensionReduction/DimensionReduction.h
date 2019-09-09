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

#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class DimensionReduction : public Debug {

  public:
    DimensionReduction();
    ~DimensionReduction();

    inline int setSEParameters(std::string &Affinity,
                               float Gamma,
                               std::string &EigenSolver) {
      se_Affinity = Affinity;
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
                                std::string &Dissimilarity) {
      mds_Metric = Metric;
      mds_Init = Init;
      mds_MaxIteration = MaxIteration;
      mds_Verbose = Verbose;
      mds_Epsilon = Epsilon;
      mds_Dissimilarity = Dissimilarity;
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
      modulePath_ = modulePath;
      return 0;
    }

    inline int setInputModuleName(const std::string &moduleName) {
      moduleName_ = moduleName;
      return 0;
    }

    inline int setInputFunctionName(const std::string &functionName) {
      functionName_ = functionName;
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

    std::string modulePath_;
    std::string moduleName_;
    std::string functionName_;
    SimplexId numberOfRows_;
    SimplexId numberOfColumns_;
    int method_;
    int numberOfComponents_;
    int numberOfNeighbors_;
    int randomState_;
    void *matrix_;
    std::vector<std::vector<double>> *embedding_;
    char majorVersion_;
  };
} // namespace ttk
