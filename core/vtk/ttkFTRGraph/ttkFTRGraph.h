/// \sa ttk::ftr::FTMTree

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

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);

  vtkSetMacro(OffsetFieldId, int);
  vtkGetMacro(OffsetFieldId, int);

  void SetSingleSweep(const bool ss) {
    params_.singleSweep = ss;
    Modified();
  }

  bool GetSingleSweep(void) const {
    return params_.singleSweep;
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

  void SetSampling(int lvl) {
    params_.samplingLvl = lvl;
    Modified();
  }

  int GetSuperArcSamplingLevel(void) const {
    return params_.samplingLvl;
  }

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
  std::string ScalarField{};
  bool ForceInputOffsetScalarField{};
  std::string InputOffsetScalarFieldName{};
  int ScalarFieldId{};
  int OffsetFieldId{};

  ttk::ftr::Params params_{};

  vtkDataSet *mesh_{};
  ttk::Triangulation *triangulation_{};
  vtkDataArray *inputScalars_{};
  std::vector<ttk::ftr::idVertex> offsets_{};
};
