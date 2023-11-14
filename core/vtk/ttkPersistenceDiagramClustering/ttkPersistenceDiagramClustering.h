/// \ingroup vtk
/// \class ttkPersistenceDiagramClustering
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \brief TTK processing package for the computation of Wasserstein barycenters
/// and K-Means clusterings of a set of persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramClustering/">Persistence
///   Diagram Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramDistance/">Persistence
///   Diagram Distance example</a> \n
///

#pragma once

// VTK includes
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkNew.h>

// VTK Module
#include <ttkPersistenceDiagramClusteringModule.h>

// TTK includes
#include <PersistenceDiagramClustering.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class vtkUnstructuredGrid;

class TTKPERSISTENCEDIAGRAMCLUSTERING_EXPORT ttkPersistenceDiagramClustering
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramClustering {

public:
  static ttkPersistenceDiagramClustering *New();

  vtkTypeMacro(ttkPersistenceDiagramClustering, ttkAlgorithm);

  enum class DISPLAY {
    COMPACT = 0,
    STARS = 1,
    MATCHINGS = 2,
  };

  enum class METHOD {
    PROGRESSIVE = 0,
    AUCTION = 1,
  };

  vtkSetMacro(WassersteinMetric, int);
  vtkGetMacro(WassersteinMetric, int);

  vtkSetMacro(UseProgressive, bool);
  vtkGetMacro(UseProgressive, bool);

  vtkSetMacro(TimeLimit, double);
  vtkGetMacro(TimeLimit, double);

  void SetAlpha(const double alpha) {
    this->Alpha = std::min(std::abs(alpha), 1.0);
    Modified();
  }
  void SetAntiAlpha(const double antiAlpha) {
    SetAlpha(1.0 - antiAlpha);
  }
  vtkGetMacro(Alpha, double);

  vtkSetMacro(DeltaLim, double);
  vtkGetMacro(DeltaLim, double);

  vtkSetMacro(Lambda, double);
  vtkGetMacro(Lambda, double);

  vtkSetMacro(NumberOfClusters, int);
  vtkGetMacro(NumberOfClusters, int);

  vtkSetMacro(UseAccelerated, bool);
  vtkGetMacro(UseAccelerated, bool);

  vtkSetMacro(UseKmeansppInit, bool);
  vtkGetMacro(UseKmeansppInit, bool);

  vtkSetMacro(ForceUseOfAlgorithm, bool);
  vtkGetMacro(ForceUseOfAlgorithm, bool);

  vtkSetMacro(Deterministic, bool);
  vtkGetMacro(Deterministic, bool);

  vtkSetMacro(PairTypeClustering, int);
  vtkGetMacro(PairTypeClustering, int);

  void SetSpacing(double spacing) {
    this->Spacing = spacing;
    Modified();
    if(!intermediateDiagrams_.empty()) {
      // skip clustering computation only if done at least once before
      this->needUpdate_ = false;
    }
  }
  vtkGetMacro(Spacing, double);

  void SetDisplayMethod(int displayMethod) {
    this->DisplayMethod = static_cast<DISPLAY>(displayMethod);
    Modified();
    if(!intermediateDiagrams_.empty()) {
      // skip clustering computation only if done at least once before
      this->needUpdate_ = false;
    }
  }
  vtkGetEnumMacro(DisplayMethod, DISPLAY);

  vtkSetMacro(UseAdditionalPrecision, bool);
  vtkGetMacro(UseAdditionalPrecision, bool);

  vtkSetMacro(DistanceWritingOptions, int);
  vtkGetMacro(DistanceWritingOptions, int);

  vtkSetMacro(UseInterruptible, bool);
  vtkGetMacro(UseInterruptible, bool);

  ttkSetEnumMacro(Method, METHOD);
  vtkGetEnumMacro(Method, METHOD);

  vtkSetMacro(NonMatchingWeight, double);
  vtkGetMacro(NonMatchingWeight, double);

protected:
  ttkPersistenceDiagramClustering();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  void Modified() override;

  void outputClusteredDiagrams(
    vtkMultiBlockDataSet *output,
    const std::vector<vtkUnstructuredGrid *> &diagsVTU,
    const std::vector<ttk::DiagramType> &diags,
    const std::vector<std::vector<std::vector<ttk::MatchingType>>>
      &matchingsPerCluster,
    const std::vector<int> &inv_clustering,
    const DISPLAY dm,
    const double spacing,
    const double max_persistence) const;
  void outputCentroids(
    vtkMultiBlockDataSet *output,
    const std::vector<ttk::DiagramType> &final_centroids,
    const std::vector<std::vector<std::vector<ttk::MatchingType>>>
      &matchingsPerCluster,
    vtkUnstructuredGrid *const someInputDiag,
    const DISPLAY dm,
    const double spacing,
    const double max_persistence) const;
  void outputMatchings(
    vtkMultiBlockDataSet *output,
    const size_t nClusters,
    const std::vector<ttk::DiagramType> &diags,
    const std::vector<std::vector<std::vector<ttk::MatchingType>>>
      &matchingsPerCluster,
    const std::vector<ttk::DiagramType> &centroids,
    const std::vector<int> &inv_clustering,
    const ttkPersistenceDiagramClustering::DISPLAY dm,
    const double spacing,
    const double max_persistence) const;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::vector<ttk::DiagramType> intermediateDiagrams_{};
  std::vector<std::vector<std::vector<ttk::MatchingType>>> all_matchings_{};
  std::vector<ttk::DiagramType> final_centroids_{};
  std::vector<int> inv_clustering_{};

  double Spacing{1.0};
  double max_dimension_total_{};

  DISPLAY DisplayMethod{DISPLAY::COMPACT};
  METHOD Method{METHOD::PROGRESSIVE};
  bool needUpdate_{true};
};
