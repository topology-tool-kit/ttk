/// \sa ttk::ftm::FTMTree
#ifndef _VTK_FTMTREE__H
#define _VTK_FTMTREE__H

// VTK includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
// Data array
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>

// Unused ? (compile without these on my computer)
// #include <vtkFiltersCoreModule.h>
// #include <vtkInformation.h>
// #include <vtkInformationVector.h>
// #include <vtkLine.h>
// #include <vtkObjectFactory.h>
// #include <vtkType.h>
// #include <vtkUnstructuredGrid.h>

// VTK module
#include <ttkFTMTreeModule.h>

// ttk code includes
#include <FTMTree.h>
#include <ttkAlgorithm.h>
#include <ttkFTMStructures.h>

class TTKFTMTREE_EXPORT ttkFTMTree : public ttkAlgorithm,
                                     protected ttk::ftm::FTMTree {

public:
  static ttkFTMTree *New();

  vtkTypeMacro(ttkFTMTree, ttkAlgorithm);

  // end of default ttk setters

  // used by command line
  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);
  vtkSetMacro(OffsetFieldId, int);
  vtkGetMacro(OffsetFieldId, int);

  // Parameters uses a structure, we can't use vtkMacro on them
  void SetTreeType(const int type) {
    params_.treeType = (ttk::ftm::TreeType)type;
    Modified();
  }

  ttk::ftm::TreeType GetTreeType(void) const {
    return params_.treeType;
  }

  void SetWithSegmentation(const bool segm) {
    params_.segm = segm;
    Modified();
  }

  bool GetWithSegmentation(void) const {
    return params_.segm;
  }

  void SetWithNormalize(const bool norm) {
    params_.normalize = norm;
    Modified();
  }

  bool GetWithNormalize(void) const {
    return params_.normalize;
  }

  void SetWithAdvStats(const bool adv) {
    params_.advStats = adv;
    Modified();
  }

  bool GetWithAdvStats(void) const {
    return params_.advStats;
  }

  void SetSuperArcSamplingLevel(int lvl) {
    params_.samplingLvl = lvl;
    Modified();
  }

  int GetSuperArcSamplingLevel(void) const {
    return params_.samplingLvl;
  }

  int setupTriangulation();
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

  // vtkDataSetAlgorithm methods

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

protected:
  ttkFTMTree();
  ~ttkFTMTree() override = default;

  void identify(vtkDataSet *ds) const;

private:
  bool ForceInputOffsetScalarField = false;
  std::string InputOffsetScalarFieldName;
  std::string ScalarFieldName;
  int ScalarFieldId = -1;
  int OffsetFieldId = -1;

  ttk::ftm::Params params_;

  int nbCC_;
  std::vector<vtkSmartPointer<vtkDataSet>> connected_components_;
  std::vector<ttk::Triangulation *> triangulation_;
  std::vector<ttk::ftm::LocalFTM> ftmTree_;
  std::vector<vtkDataArray *> inputScalars_;
  std::vector<std::vector<ttk::SimplexId>> offsets_;

  bool hasUpdatedMesh_;
};

#endif // _VTK_FTMTREE__H
