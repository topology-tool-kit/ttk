/// \ingroup vtk
/// \class ttkFTMTree
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date June 2017.
///
/// \sa ttk::ftm::FTMTree
///
/// \brief TTK VTK-filter for the computation of merge and contour trees.
///
/// The computation of the Merge / Contour tree done by this package is done in
/// parallel if TTK_ENABLE_OPENMP is set to ON, using a task based approch
/// described in the article mention below.
/// The VTK wrapper will first call a connectivity filter, and then call
/// a contour / merge tree computation for each connected components. The final
/// tree is then aggregated.
///
/// \param Input Input scalar field, either 2D or 3D, regular
/// grid or triangulation (vtkDataSet)
/// \param TreeType the Type of three to Compute:\n
/// * Join Tree (leaves corresponds to minima of the scalar field)
/// * Split Tree (leaves corresponds to maxima of the scalar field)
/// * Contour Tree (combination of both)
/// * JoinSplit (compute both merge trees but do not combine them (advanced))
/// \param Segmentation control wethear or not the output should be augmented
/// with the segmentation.
/// \param SuperArcSamplingLevel control the number of subdivision of each
/// superarc. Intermediate point will be located on the barycenter of the
/// corresponding portion of vertex.
/// \param Output the output of this filter is composed of:\n
/// 1. The nodes of the tree
/// 2. The arcs of the tree
/// 3. The semgentation of the initial dataset
/// The structure of the tree (Nodes+Arcs) have a concept of nodeId, wich is
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
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
///   example</a>
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeTemporalReduction/">Merge
///   Tree Temporal Reduction</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n

#pragma once

// VTK includes
#include <vtkSmartPointer.h>

// VTK module
#include <ttkFTMTreeModule.h>

// ttk code includes
#include <FTMTree.h>
#include <ttkAlgorithm.h>
#include <ttkFTMStructures.h>

class vtkDataSet;

class TTKFTMTREE_EXPORT ttkFTMTree : public ttkAlgorithm {

public:
  static ttkFTMTree *New();

  vtkTypeMacro(ttkFTMTree, ttkAlgorithm);

  /// @brief the offset array to use for simulation of simplicity
  /// @{
  vtkGetMacro(ForceInputOffsetScalarField, bool);
  vtkSetMacro(ForceInputOffsetScalarField, bool);
  /// @}

  // Parameters uses a structure, we can't use vtkMacro on them

  /// @brief the type of tree to compute (Join, Split, Contour, JoinSplit)
  /// @{
  void SetTreeType(const int type) {
    params_.treeType = (ttk::ftm::TreeType)type;
    Modified();
  }
  ttk::ftm::TreeType GetTreeType(void) const {
    return params_.treeType;
  }
  /// @}

  /// @brief control if the output should contains the segmentation information
  /// @{
  void SetWithSegmentation(const bool segm) {
    params_.segm = segm;
    Modified();
  }
  bool GetWithSegmentation(void) const {
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
  bool GetWithNormalize(void) const {
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
  bool GetWithAdvStats(void) const {
    return params_.advStats;
  }
  /// @}

  /// @brief control the sampling level of the superarc. By default: 0
  /// @{
  void SetSuperArcSamplingLevel(int lvl) {
    params_.samplingLvl = lvl;
    Modified();
  }
  int GetSuperArcSamplingLevel(void) const {
    return params_.samplingLvl;
  }
  /// @}

  int preconditionTriangulation();
  int getScalars();
  int getOffsets();

  int getSkeletonNodes(vtkUnstructuredGrid *outputSkeletonNodes);

  int addDirectSkeletonArc(const ttk::ftm::idSuperArc arcId,
                           const int cc,
                           vtkPoints *points,
                           vtkUnstructuredGrid *skeletonArcs,
                           ttk::ftm::ArcData &arcData);

  int addSampledSkeletonArc(const ttk::ftm::idSuperArc arcId,
                            const int cc,
                            vtkPoints *points,
                            vtkUnstructuredGrid *skeletonArcs,
                            ttk::ftm::ArcData &arcData);

  int addCompleteSkeletonArc(const ttk::ftm::idSuperArc arcId,
                             const int cc,
                             vtkPoints *points,
                             vtkUnstructuredGrid *skeletonArcs,
                             ttk::ftm::ArcData &arcData);

  int getSkeletonArcs(vtkUnstructuredGrid *outputSkeletonArcs);

  int getSegmentation(vtkDataSet *outputSegmentation);

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
  void printCSVStats();
  void printCSVTree(const ttk::ftm::FTMTree_MT *const tree) const;
#endif

protected:
  ttkFTMTree();

  // vtkDataSetAlgorithm methods
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  void identify(vtkDataSet *ds) const;

private:
  bool ForceInputOffsetScalarField = false;
  ttk::ftm::Params params_;

  int nbCC_;
  std::vector<vtkSmartPointer<vtkDataSet>> connected_components_;
  std::vector<ttk::Triangulation *> triangulation_;
  std::vector<ttk::ftm::LocalFTM> ftmTree_;
  std::vector<vtkDataArray *> inputScalars_;
  std::vector<std::vector<ttk::SimplexId>> offsets_;
};
