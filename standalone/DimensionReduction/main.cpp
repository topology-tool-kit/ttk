/// \author Alexandre Talon <alexandre.talon@lip6.fr>.
/// \date November 2023
///
/// \brief dummy program example.

// TTK Includes
#include <CommandLineParser.h>
#include <ttkDimensionReduction.h>

// VTK Includes
#include <vtkCellData.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

int main(int argc, char **argv) {
  std::vector<std::string> inputFilePaths;
  std::vector<std::string> inputArrayNames;
  std::string outputArrayName{"2d-projection"};
  std::string outputPathPrefix{"output"};
  bool listArrays{false};

  // Variables for the module options.
  int backend{0};
  bool isInputDistMat{false};
  int nbComp{2}, nbNeighb{5};
  bool isDeterministic{true};

  // Spectral Embedding
  std::string se_Affinity{"nearest_neighbors"};
  double se_Gamma{1};
  std::string se_EigenSolver{"None"};

  // Locally Linear Embedding
  double lle_Regularization{1e-3};
  std::string lle_EigenSolver{"auto"};
  double lle_Tolerance{1e-3};
  int lle_MaxIteration{300};
  std::string lle_Method{"standard"};
  double lle_HessianTolerance{1e-3};
  double lle_ModifiedTolerance{1e-3};
  std::string lle_NeighborsAlgorithm{"auto"};

  // Multi-Dimensional Scaling
  bool mds_Metric{true};
  int mds_Init{4};
  int mds_MaxIteration{300};
  int mds_Verbose{0};
  double mds_Epsilon{0};

  // t-distributed Stochastic Neighbor Embedding
  double tsne_Perplexity{30};
  double tsne_Exaggeration{12};
  double tsne_LearningRate{200};
  int tsne_MaxIteration{1000};
  int tsne_MaxIterationProgress{300};
  double tsne_GradientThreshold{1e-7};
  std::string tsne_Metric{"euclidean"};
  std::string tsne_Init{"pca"};
  int tsne_Verbose{0};
  std::string tsne_Method{"barnes_hut"};
  double tsne_Angle{0.5};

  // Isomap Embedding
  std::string iso_EigenSolver{"auto"};
  float iso_Tolerance{1e-3};
  int iso_MaxIteration{300};
  std::string iso_PathMethod{"auto"};
  std::string iso_NeighborsAlgorithm{"auto"};
  std::string iso_Metric{"euclidean"};

  // Principal Component Analysis
  bool pca_Copy{true};
  bool pca_Whiten{false};
  std::string pca_SVDSolver{"auto"};
  double pca_Tolerance{0};
  std::string pca_MaxIteration{"auto"};

  // TopoMap
  int topomap_AngularSampleNb{2};
  int topomap_CheckMST{false};
  int topomap_Strategy{0};

  // Set program variables based on command line arguments
  {
    ttk::CommandLineParser parser;

    // Standard options and arguments
    parser.setArgument("i", &inputFilePaths, "Input data-sets (*.csv)", false);
    parser.setArgument("a", &inputArrayNames, "Input array names", true);
    parser.setArgument(
      "o", &outputPathPrefix, "Output file prefix (*.csv)", true);
    parser.setOption("l", &listArrays, "List available arrays");

    // Declare custom arguments and options
    parser.setArgument(
      "B", &backend,
      "Method (0: Spectral Embedding, 1: Locally Linear Embedding, 2. "
      "Multi-Dimensional Scaling, 3: t-distributed Stochastic Neighbor "
      "Embedding, 4: Isomap Embedding, 5: Principal Component Analysis, 6: "
      "TopoMap (IEEE VIS 2020))",
      true);
    parser.setArgument("Dist", &isInputDistMat,
                       "1 if the input is a distance matrix, 0 otherwise.",
                       true);
    parser.setArgument("C", &nbComp, "Number of components", true);
    parser.setArgument("N", &nbNeighb, "Number of neighbors", true);
    parser.setArgument("Det", &isDeterministic,
                       "1 if we want a deterministic output, 0 otherwise.",
                       true);

    parser.setArgument("O", &outputArrayName, "Output array name", true);

    // Spectral Embedding
    parser.setArgument(
      "se_affinity", &se_Affinity,
      "SE: Specifies how to construct the affinity matrix if the input is not "
      "a distance matrix [nearest_neighbors, rbf].",
      true);
    parser.setArgument("se_gamma", &se_Gamma,
                       "SE: Kernel coefficient if using rbf for affinity.",
                       true);
    parser.setArgument("se_solver", &se_EigenSolver,
                       "SE: The type of eigen solver [arpack, lobpcg, amg]",
                       true);

    // Locally Linear Embedding
    parser.setArgument(
      "lle_regul", &lle_Regularization, "Lle: Regularization value", true);
    parser.setArgument("lle_solver", &lle_EigenSolver,
                       "Lle: The type of eigen solver [auto, arpack, dense]",
                       true);
    parser.setArgument("lle_tolerance", &lle_Tolerance,
                       "Lle: The tolerance value (only for arpack solver).",
                       true);
    parser.setArgument("lle_iterThreshold", &lle_MaxIteration,
                       "Lle: The iteration threshold.", true);
    parser.setArgument(
      "lle_method", &lle_Method,
      "Isomap: The embedding method [standard, hessian, modified, ltsa].",
      true);
    parser.setArgument("lle_hessianTolerance", &lle_HessianTolerance,
                       "Lle: The tolerance value for the Hessian method.",
                       true);
    parser.setArgument("lle_modifiedTolerance", &lle_ModifiedTolerance,
                       "Lle: The tolerance value for the Modified method.",
                       true);
    parser.setArgument(
      "lle_nnAlgo", &lle_NeighborsAlgorithm,
      "Lle: The neighbors algorithm [auto, brute, kd_tree, ball_tree].", true);

    // Multi-Dimensional Scaling
    parser.setArgument("mds_metric", &mds_Metric,
                       "MDS: To perform metric or nonmetric MDS.", true);
    parser.setArgument(
      "mds_init", &mds_Init,
      "MDS: Number of times the SMACOF algorithm will be run with different "
      "initializations. The final results will be the best output of the runs, "
      "determined by the run with the smallest final stress.",
      true);
    parser.setArgument("mds_maxIter", &mds_MaxIteration,
                       "MDS: Maximum number of iterations of the SMACOF "
                       "algorithm for a single run.",
                       true);
    parser.setArgument(
      "mds_verb", &mds_Verbose, "MDS: Level of verbosity.", true);
    parser.setArgument("mds_epsilon", &mds_Epsilon,
                       "MDS: Relative tolerance with respect to stress at "
                       "which to declare convergence.",
                       true);

    // t-distributed Stochastic Neighbor Embedding
    parser.setArgument(
      "tsne_perplexity", &tsne_Perplexity,
      "t-sne: The perplexity is related to the number of nearest neighbors "
      "that is used in other manifold learning algorithms. Larger datasets "
      "usually require a larger perplexity. Consider selecting a value between "
      "5 and 50.",
      true);
    parser.setArgument(
      "tsne_exaggeration", &tsne_Exaggeration,
      "t-sne: The early exxageration controls how tight natural clusters in "
      "the original space are in the embedded space and how much space will be "
      "between them. For larger values, the space between natural clusters "
      "will be larger in the embedded space.",
      true);
    parser.setArgument(
      "tsne_learningRate", &tsne_LearningRate,
      "t-sne: If the learning rate is too high, the data may look like a "
      "‘ball’ with any point approximately equidistant from its nearest "
      "neighbours. If the learning rate is too low, most points may look "
      "compressed in a dense cloud with few outliers. If the cost function "
      "gets stuck in a bad local minimum increasing the learning rate may "
      "help. Usually: [10.0, 1000.0]",
      true); // TODO allow 'auto'?
    parser.setArgument("tsne_maxIter", &tsne_MaxIteration,
                       "t-sne: Maximum number of iterations for the "
                       "optimization. Should be at least 250.",
                       true);
    parser.setArgument(
      "tsne_maxIterNoProgress", &tsne_MaxIterationProgress,
      "t-sne: Maximum number of iterations without progress before we abort "
      "the optimization, used after 250 initial iterations with early "
      "exaggeration. This value is rounded to the next multiple of 50.",
      true);
    parser.setArgument("tsne_gradientThresh", &tsne_GradientThreshold,
                       "t-sne: If the gradient norm is below this threshold, "
                       "the optimization will be stopped.",
                       true);
    parser.setArgument(
      "tsne_metric", &tsne_Metric,
      "t-sne: The metric to use when calculating distance between instances in "
      "a feature array. [euclidean, ... (see "
      "https://scikit-learn.org/stable/modules/generated/"
      "sklearn.metrics.pairwise_distances.html#sklearn.metrics.pairwise_"
      "distances)].",
      true);
    parser.setArgument(
      "tsne_init", &tsne_Init,
      "t-sne: Initialization of embedding. PCA initialization cannot be used "
      "with precomputed distances and is usually more globally stable than "
      "random initialization.",
      true);
    parser.setArgument(
      "tsne_verb", &tsne_Verbose, "t-sne: Verbosity level.", true);
    parser.setArgument(
      "tsne_method", &tsne_Method,
      "t-sne: By default the gradient calculation algorithm uses Barnes-Hut "
      "approximation running in O(NlogN) time. exact will run in time O(N^2) "
      "time. The exact algorithm should be used when nearest-neighbor errors "
      "need to be better than 3%. [barnes_hut, exact]",
      true);
    parser.setArgument(
      "tsne_angle", &tsne_Angle,
      "t-sne: Only used if method=’barnes_hut’. This is the trade-off between "
      "speed and accuracy for Barnes-Hut T-SNE. For more details, see "
      "https://scikit-learn.org/stable/modules/generated/"
      "sklearn.manifold.TSNE.html#sklearn.manifold.TSNE. [",
      true);

    // Isomap Embedding
    parser.setArgument("iso_solver", &iso_EigenSolver,
                       "Isomap: The type of eigen solver [auto, arpack, dense]",
                       true);
    parser.setArgument("iso_tolerance", &iso_EigenSolver,
                       "Isomap: The tolerance value (only for arpack solver).",
                       true);
    parser.setArgument(
      "iso_maxIter", &iso_MaxIteration,
      "Isomap: The maximum number of iteration (only for arpack solver).",
      true);
    parser.setArgument("iso_pathMethod", &iso_PathMethod,
                       "Isomap: The shortest path method [auto, FD, D] "
                       "(Floyd-Warshall / Dijkstra).",
                       true);
    parser.setArgument(
      "iso_nnAlgo", &iso_NeighborsAlgorithm,
      "Isomap: The neighbors algorithm [auto, brute, kd_tree, ball_tree].",
      true);
    parser.setArgument("iso_metric", &iso_Metric,
                       "Isomap: The metric used [euclidean, ... (see "
                       "https://scikit-learn.org/stable/modules/generated/"
                       "sklearn.metrics.pairwise_distances.html#sklearn."
                       "metrics.pairwise_distances)].",
                       true);

    // Principal Component Analysis
    parser.setArgument(
      "pca_copy", &pca_Copy,
      "PCA: Does not overwrite the data passed to fit if the value is true.",
      true);
    parser.setArgument("pca_whiten", &pca_Whiten,
                       "PCA: Whiten option (see "
                       "https://scikit-learn.org/stable/modules/generated/"
                       "sklearn.decomposition.PCA.html for details).",
                       true);
    parser.setArgument("pca_solver", &pca_SVDSolver,
                       "PCA: SVD Solver [auto, full, arpack, randomized].",
                       true);
    parser.setArgument("pca_tolerance", &pca_Tolerance,
                       "PCA: The tolerance value (only for arpack solver).",
                       true);
    parser.setArgument(
      "pca_maxIter", &pca_Tolerance,
      "PCA: Number of iterations for the power method computed by thesolver "
      "(only for randomized solver) ['auto' or some integer].",
      true);

    // TopoMap
    parser.setArgument("tm_strat", &topomap_Strategy,
                       "TopoMap: Strategy for TopoMap (0: Kruskal, 1: Prim)",
                       true);

    parser.setArgument("tm_angles", &topomap_AngularSampleNb,
                       "TopoMap: Number of angular samples (topomp)", true);
    parser.setArgument(
      "tm_check", &topomap_CheckMST,
      "Topomap: 1 if the MST should be checked for consistency, 0 otherwise.",
      true);

    parser.parse(argc, argv);
  }

  ttk::Debug msg;
  msg.setDebugMsgPrefix("DimensionReduction");
  auto dimRed = vtkSmartPointer<ttkDimensionReduction>::New();

  // Pass custom arguments and options to the module
  dimRed->SetSelectFieldsWithRegexp(true);
  dimRed->SetRegexpString(".*");

  dimRed->SetInputIsADistanceMatrix(isInputDistMat);
  dimRed->SetMethod(backend);
  dimRed->SetNumberOfComponents(nbComp);
  dimRed->SetNumberOfNeighbors(nbNeighb);
  dimRed->SetIsDeterministic(isDeterministic);

  dimRed->Setse_Affinity(se_Affinity);
  dimRed->Setse_Gamma(se_Gamma);
  dimRed->Setse_EigenSolver(se_EigenSolver);

  dimRed->Setlle_Regularization(lle_Regularization);
  dimRed->Setlle_EigenSolver(lle_EigenSolver);
  dimRed->Setlle_Tolerance(lle_Tolerance);
  dimRed->Setlle_MaxIteration(lle_MaxIteration);
  dimRed->Setlle_Method(lle_Method);
  dimRed->Setlle_HessianTolerance(lle_HessianTolerance);
  dimRed->Setlle_ModifiedTolerance(lle_ModifiedTolerance);
  dimRed->Setlle_NeighborsAlgorithm(lle_NeighborsAlgorithm);

  dimRed->Setmds_Metric(mds_Metric);
  dimRed->Setmds_Init(mds_Init);
  dimRed->Setmds_MaxIteration(mds_MaxIteration);
  dimRed->Setmds_Verbose(mds_Verbose);
  dimRed->Setmds_Epsilon(mds_Epsilon);

  dimRed->Settsne_Perplexity(tsne_Perplexity);
  dimRed->Settsne_Exaggeration(tsne_Exaggeration);
  dimRed->Settsne_LearningRate(tsne_LearningRate);
  dimRed->Settsne_MaxIteration(tsne_MaxIteration);
  dimRed->Settsne_MaxIterationProgress(tsne_MaxIterationProgress);
  dimRed->Settsne_GradientThreshold(tsne_GradientThreshold);
  dimRed->Settsne_Metric(tsne_Metric);
  dimRed->Settsne_Init(tsne_Init);
  dimRed->Settsne_Verbose(tsne_Verbose);
  dimRed->Settsne_Method(tsne_Method);
  dimRed->Settsne_Angle(tsne_Angle);

  dimRed->Setiso_EigenSolver(iso_EigenSolver);
  dimRed->Setiso_Tolerance(iso_Tolerance);
  dimRed->Setiso_MaxIteration(iso_MaxIteration);
  dimRed->Setiso_PathMethod(iso_PathMethod);
  dimRed->Setiso_NeighborsAlgorithm(iso_NeighborsAlgorithm);
  dimRed->Setiso_Metric(iso_Metric);

  dimRed->Setpca_Copy(pca_Copy);
  dimRed->Setpca_Whiten(pca_Whiten);
  dimRed->Setpca_SVDSolver(pca_SVDSolver);
  dimRed->Setpca_Tolerance(pca_Tolerance);
  dimRed->Setpca_MaxIteration(pca_MaxIteration);

  dimRed->Settopomap_Strategy(topomap_Strategy);
  dimRed->Settopomap_AngularSampleNb(topomap_AngularSampleNb);
  dimRed->Settopomap_CheckMST(topomap_CheckMST);

  // ---------------------------------------------------------------------------
  // Read input vtkDataObjects (optionally: print available arrays)
  // ---------------------------------------------------------------------------
  for(size_t i = 0; i < inputFilePaths.size(); i++) {
    auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetFileName(inputFilePaths[i].data());
    reader->DetectNumericColumnsOn();
    reader->Update();
    // check if input vtkDataObject was successfully read
    auto inputDataObject = reader->GetOutput();

    if(!inputDataObject) {
      msg.printErr("Unable to read input file `" + inputFilePaths[i] + "' :(");
      return 1;
    }

    auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject);

    // if requested print list of arrays, otherwise proceed with execution
    if(listArrays) {
      msg.printMsg(inputFilePaths[i] + ":");
      if(inputAsVtkDataSet) {
        // Point Data
        msg.printMsg("  PointData:");
        auto pointData = inputAsVtkDataSet->GetPointData();
        for(int j = 0; j < pointData->GetNumberOfArrays(); j++)
          msg.printMsg("    - " + std::string(pointData->GetArrayName(j)));

        // Cell Data
        msg.printMsg("  CellData:");
        auto cellData = inputAsVtkDataSet->GetCellData();
        for(int j = 0; j < cellData->GetNumberOfArrays(); j++)
          msg.printMsg("    - " + std::string(cellData->GetArrayName(j)));
      } else {
        msg.printErr("Unable to list arrays on file `" + inputFilePaths[i]
                     + "'");
        return 1;
      }
    } else {
      // feed input object to ttkHelloWorld filter
      dimRed->SetInputDataObject(0, reader->GetOutput());
    }
  }

  // terminate program if it was just asked to list arrays
  if(listArrays) {
    return 0;
  }

  dimRed->Update();

  // If output prefix is specified then write all output objects to disk
  if(!outputPathPrefix.empty()) {
    for(int i = 0; i < dimRed->GetNumberOfOutputPorts(); i++) {
      auto output = dimRed->GetOutputPort(i);
      auto writer = vtkSmartPointer<vtkDelimitedTextWriter>::New();

      std::string outputFileName
        = outputPathPrefix + "_port_" + std::to_string(i) + ".csv";
      msg.printMsg("Writing output file `" + outputFileName + "'...");
      writer->SetInputConnection(output);
      writer->SetFileName(outputFileName.data());
      writer->Update();
    }
  }

  return 0;
}
