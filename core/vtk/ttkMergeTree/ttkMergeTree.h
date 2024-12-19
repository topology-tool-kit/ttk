/// \ingroup vtk
/// \class ttkMergeTree
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date June 2017.
///
/// \sa ttk::ftm::FTMTree
///
/// \brief TTK filter for the computation of merge trees.
///
/// The computation of the merge tree done by this package is done in
/// parallel if TTK_ENABLE_OPENMP is set to ON, using a task based approach
/// described in the article mention below.
/// The VTK wrapper will first call a connectivity filter, and then call
/// a merge tree computation for each connected components. The final
/// tree is then aggregated.
///
/// \param Input Input scalar field, either 2D or 3D, regular
/// grid or triangulation (vtkDataSet)
/// \param TreeType the Type of three to Compute:\n
/// * Join Tree (leaves corresponds to minima of the scalar field)
/// * Split Tree (leaves corresponds to maxima of the scalar field)
/// \param Segmentation control wethear or not the output should be augmented
/// with the segmentation.
/// \param SuperArcSamplingLevel control the number of subdivision of each
/// superarc. Intermediate point will be located on the barycenter of the
/// corresponding portion of vertex.
/// \param Output the output of this filter is composed of:\n
/// 1. The nodes of the tree
/// 2. The arcs of the tree
/// 3. The semgentation of the initial dataset
/// The structure of the tree (Nodes+Arcs) have a concept of nodeId, which is
/// an id that is consistent between execution if SetWithNormalize is set to
/// True. The downNodeId of an arc is its starting node (directed towards the
/// leaves as the computation starts here) and the upNodeId it the ending node,
/// in direction of the Root of the tree.
/// The segmentation also contains some basics metrics like the size of each
/// region (RegionSpan) or its number of vertex (RegionSize)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \b Related \b publication \n
/// "Task-based Augmented Merge Trees with Fibonacci Heaps" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// 2017 IEEE 7th Symposium on Large Data Analysis and Visualization (LDAV),
/// doi: 10.1109/LDAV.2017.8231846. \n
/// "Task-based augmented contour trees with fibonacci heaps" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// IEEE Transactions on Parallel and Distributed Systems, Volume 30, Issue 8,
/// Pages 1889-1905
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/ctBones/">CT Bones
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeFeatureTracking/">Merge
///   Tree Feature Tracking example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeTemporalReduction/">Merge
///   Tree Temporal Reduction</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeWAE/">Merge
///   tree Wasserstein Auto-Encoder example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeExTreeM/">Merge
///   Trees Via ExTreeM</a> \n

#pragma once

// VTK includes
#include <vtkSmartPointer.h>

// VTK module
#include <ttkMergeTreeModule.h>

// ttk code includes
#include <PathCompression.h>
#include <ttkAlgorithm.h>
#include <ttkMergeTreeBase.h>
#include <ttkMergeTreeStructures.h>

class vtkDataSet;

class TTKMERGETREE_EXPORT ttkMergeTree : public ttkAlgorithm,
                                         public ttkMergeTreeBase {

public:
  enum class BACKEND {
    FTM = 0,
    EXTREEM = 1,
  };

  static ttkMergeTree *New();

  vtkTypeMacro(ttkMergeTree, ttkAlgorithm);

  /// @brief the offset array to use for simulation of simplicity
  /// @{
  vtkGetMacro(ForceInputOffsetScalarField, bool);
  vtkSetMacro(ForceInputOffsetScalarField, bool);
  /// @}

  /// @brief the backend to use for computations.
  /// @{
  vtkGetMacro(Backend, int);
  vtkSetMacro(Backend, int);
  /// @}

  // Parameters uses a structure, we can't use vtkMacro on them

  /// @brief the type of tree to compute (Join, Split, Contour, JoinSplit)
  /// @{
  void SetTreeType(const int type) {
    params_.treeType = (ttk::ftm::TreeType)type;
    Modified();
  }
  /// @}

  /// @brief control if the output should contains the segmentation information
  /// @{
  void SetWithSegmentation(const bool segm) {
    params_.segm = segm;
    Modified();
  }
  bool GetWithSegmentation() const {
    return params_.segm;
  }
  /// @}

  /// @brief if true, a post process pass will ensure NodesId have a
  /// deterministic order
  /// @{
  void SetWithNormalize(const bool norm) {
    params_.normalize = norm;
    Modified();
  }
  bool GetWithNormalize() const {
    return params_.normalize;
  }
  /// @}

  /// @brief Compute additional information on the segmentation
  /// like the span and size (in nb of vertex) of each region
  /// @{
  void SetWithAdvStats(const bool adv) {
    params_.advStats = adv;
    Modified();
  }
  bool GetWithAdvStats() const {
    return params_.advStats;
  }
  /// @}

  /// @brief control the sampling level of the superarc. By default: 0
  /// @{
  void SetSuperArcSamplingLevel(int lvl) {
    params_.samplingLvl = lvl;
    Modified();
  }
  int GetSuperArcSamplingLevel() const {
    return params_.samplingLvl;
  }
  /// @}

protected:
  ttkMergeTree();

  int getOffsets();
  int getScalars();

  int preconditionTriangulation();

  // vtkDataSetAlgorithm methods
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int Backend{(int)BACKEND::FTM};
};
