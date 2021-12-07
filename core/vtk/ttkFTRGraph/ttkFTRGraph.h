/// \ingroup vtk
/// \class ttkFTRGraph
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date June 2017.
///
/// \sa ttk::ftm::FTRGraph
///
/// \brief TTK VTK-filter for the computation of Reeb Graphs
///
/// The computation of the Reeb graph done by this package is done in
/// parallel if TTK_ENABLE_OPENMP is set to ON, using a task based approch
/// described in the article mention below.
///
/// \param Input Input scalar field, either 2D or 3D, regular
/// grid or triangulation (vtkDataSet)
/// \param SingleSweep control if the computation should start from both minima
/// and maxima simultaneously. If you encouter troubled with FTR, you should try
/// to use the single sweep. It is slower but may be more robust.
/// \param Segmentation control wethear or not the output should be augmented
/// with the segmentation.
/// \param SuperArcSamplingLevel control the number of subdivision
/// of each superarc. Intermediate point will be located on the barycenter of
/// the corresponding portion of vertex.
/// \param Output the output of this filter
/// is composed of:\n
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
/// "Task-based Augmented Reeb Graphs with Dynamic ST-Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// EGPGV19: Eurographics Symposium on Parallel Graphics and Visualization
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n

#pragma once

// ttk code includes
#include "FTRGraph.h"
#include "Graph.h"
#include "ttkAlgorithm.h"
#include "ttkFTRGraphStructures.h"

// VTK includes
#include <vtkDataArray.h>

// VTK Module
#include <ttkFTRGraphModule.h>

class TTKFTRGRAPH_EXPORT ttkFTRGraph : public ttkAlgorithm {
public:
  static ttkFTRGraph *New();
  vtkTypeMacro(ttkFTRGraph, ttkAlgorithm);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  /// @brief control whether the computation should start from min and max
  /// (default) or use a single sweep starting only from the min (may be more
  /// robust)
  /// @{
  void SetSingleSweep(const bool ss) {
    params_.singleSweep = ss;
    Modified();
  }
  bool GetSingleSweep(void) const {
    return params_.singleSweep;
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

  /// @brief if set to true, a post processing pass will be used to enforce
  /// consistend node ids between executions
  /// @{
  void SetWithNormalize(const bool norm) {
    params_.normalize = norm;
    Modified();
  }
  bool GetWithNormalize(void) const {
    return params_.normalize;
  }
  /// @}

  /// @brief control the sampling level of the superarcs
  /// @{
  void SetSampling(int lvl) {
    params_.samplingLvl = lvl;
    Modified();
  }
  int GetSuperArcSamplingLevel(void) const {
    return params_.samplingLvl;
  }
  /// @}

  int getSkeletonNodes(const ttk::ftr::Graph &graph,
                       vtkUnstructuredGrid *outputSkeletonNodes);

  int addDirectSkeletonArc(const ttk::ftr::Graph &graph,
                           const ttk::ftr::idSuperArc arcId,
                           vtkPoints *points,
                           vtkUnstructuredGrid *skeletonArcs,
                           ttk::ftr::ArcData &arcData);

  int addSampledSkeletonArc(const ttk::ftr::Graph &graph,
                            const ttk::ftr::idSuperArc arcId,
                            vtkPoints *points,
                            vtkUnstructuredGrid *skeletonArcs,
                            ttk::ftr::ArcData &arcData);

  int addCompleteSkeletonArc(const ttk::ftr::Graph &graph,
                             const ttk::ftr::idSuperArc arcId,
                             vtkPoints *points,
                             vtkUnstructuredGrid *skeletonArcs,
                             ttk::ftr::ArcData &arcData);

  int getSkeletonArcs(const ttk::ftr::Graph &graph,
                      vtkUnstructuredGrid *outputSkeletonArcs);

  int getSegmentation(const ttk::ftr::Graph &graph,
                      vtkDataSet *outputSegmentation);

  template <typename VTK_TT, typename TTK_TT>
  int dispatch(ttk::ftr::Graph &graph);

protected:
  ttkFTRGraph();

  void identify(vtkDataSet *ds) const;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceInputOffsetScalarField{};
  ttk::ftr::Params params_{};

  vtkDataSet *mesh_{};
  ttk::Triangulation *triangulation_{};
  vtkDataArray *inputScalars_{};
  vtkDataArray *offsets_{};
};
