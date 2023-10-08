/// \ingroup vtk
/// \class ttkMergeTree
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date June 2017.
///
/// \sa ttk::ftm::FTMTree
///
/// \brief VTK-level features for merge tree representation.
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

#pragma once

// ttk code includes
#include <ttkMergeTreeStructures.h>

class ttkMergeTreeBase : virtual public ttk::Debug {

public:
  ttkMergeTreeBase();
  ~ttkMergeTreeBase() override = default;

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

  ttk::ftm::TreeType GetTreeType() const {
    return params_.treeType;
  }

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
  void printCSVStats();
  void printCSVTree(const ttk::ftm::FTMTree_MT *const tree) const;
#endif

protected:
  void identify(vtkDataSet *ds) const;

  bool ForceInputOffsetScalarField = false;
  ttk::ftm::Params params_;
  int nbCC_;
  std::vector<vtkSmartPointer<vtkDataSet>> connected_components_;
  std::vector<ttk::Triangulation *> triangulation_;
  std::vector<ttk::ftm::LocalFTM> ftmTree_;
  std::vector<vtkDataArray *> inputScalars_;
  std::vector<std::vector<ttk::SimplexId>> offsets_;
};
