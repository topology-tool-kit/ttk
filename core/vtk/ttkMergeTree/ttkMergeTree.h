/// \ingroup vtk
/// \class ttkMergeTree
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date September 2023.
///
/// \brief TTK VTK-filter for the computation of  merge trees.
///
/// \param Input Input scalar field, either 2D or 3D, regular
/// grid or triangulation (vtkDataSet)
/// \param TreeType the Type of three to Compute:\n
/// * Join Tree (leaves corresponds to minima of the scalar field)
/// * Split Tree (leaves corresponds to maxima of the scalar field)
/// \param Output the output of this filter is composed of:\n
/// 1. The nodes of the tree
/// 2. The arcs of the tree
/// 3. The semgentation of the initial dataset
///
/// This filter can be used like any other VTK filter (for instance, by using
/// the sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MergeTree/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "ExTreeM: Scalable Augmented Merge Tree Computation via Extremum Graphs"
/// \n Lukasczyk et al. \n IEEE VIS 2023
#pragma once

// VTK Module
#include <ttkMergeTreeModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <ExTreeM.h>

class vtkUnstructuredGrid;

class TTKMERGETREE_EXPORT ttkMergeTree
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ExTreeM // and we inherit from the base class
{
private:
  int Type{0};

public:
  vtkGetMacro(Type, int);
  vtkSetMacro(Type, int);

  static ttkMergeTree *New();
  vtkTypeMacro(ttkMergeTree, ttkAlgorithm);

protected:
  ttkMergeTree();
  ~ttkMergeTree() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class triangulationType = ttk::AbstractTriangulation>
  int getMergeTree(vtkUnstructuredGrid *outputSkeletonArcs,
                   std::vector<ExTreeM::Branch> &mergeTree,
                   vtkDataArray *inputScalars,
                   const triangulationType *triangulation);

  template <class triangulationType = ttk::AbstractTriangulation>
  int getMergeTreePoints(
    vtkUnstructuredGrid *outputSkeletonNodes,
    std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &persistencePairs,
    vtkDataArray *inputScalars,
    const triangulationType *triangulation);
};
