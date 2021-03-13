/// \ingroup base
/// \class ttk::ttkPersistenceDiagramBarycenter
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

#pragma once

#ifndef diagramTuple
#define diagramTuple                                                       \
  std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,            \
             ttk::CriticalType, dataType, ttk::SimplexId, dataType, float, \
             float, float, dataType, float, float, float>
#endif

#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceDiagramClusteringModule.h>

// ttk code includes
#include <PersistenceDiagramBarycenter.h>
#include <PersistenceDiagramClustering.h>

#include <ttkAlgorithm.h>

class TTKPERSISTENCEDIAGRAMCLUSTERING_EXPORT ttkPersistenceDiagramClustering
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramClustering {

public:
  static ttkPersistenceDiagramClustering *New();

  vtkTypeMacro(ttkPersistenceDiagramClustering, ttkAlgorithm);

  vtkSetMacro(WassersteinMetric, int);
  vtkGetMacro(WassersteinMetric, int);

  vtkSetMacro(UseProgressive, bool);
  vtkGetMacro(UseProgressive, bool);

  vtkSetMacro(TimeLimit, double);
  vtkGetMacro(TimeLimit, double);

  void SetAlpha(double data) {
    if(data > 0 && data <= 1) {

      Alpha = data;
    } else if(data > 1) {
      Alpha = 1;
    } else {
      Alpha = 0.001;
    }
    Modified();
  }
  vtkGetMacro(Alpha, double);

  void SetAntiAlpha(double data) {
    double alpha = 1 - data;
    SetAlpha(alpha);
  }

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
    Spacing = spacing;
    oldSpacing = spacing;
    Modified();
    if(!intermediateDiagrams_.empty()) {
      // skip clustering computation only if done at least once before
      needUpdate_ = false;
    }
  }
  vtkGetMacro(Spacing, double);

  void SetDisplayMethod(int displayMethod) {
    DisplayMethod = displayMethod;
    if(displayMethod == 0) { // compact display
      Spacing = 0;
    } else {
      Spacing = oldSpacing;
    }
    Modified();
    if(!intermediateDiagrams_.empty()) {
      // skip clustering computation only if done at least once before
      needUpdate_ = false;
    }
  }

  vtkGetMacro(DisplayMethod, bool);

  vtkSetMacro(UseAdditionalPrecision, bool);
  vtkGetMacro(UseAdditionalPrecision, bool);

  vtkSetMacro(DistanceWritingOptions, int);
  vtkGetMacro(DistanceWritingOptions, int);

  vtkSetMacro(UseInterruptible, bool);
  vtkGetMacro(UseInterruptible, bool);

  vtkSetMacro(Method, double);
  vtkGetMacro(Method, double);

protected:
  ttkPersistenceDiagramClustering();

  using diagramType = std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 double,
                                 ttk::SimplexId,
                                 double,
                                 float,
                                 float,
                                 float,
                                 double,
                                 float,
                                 float,
                                 float>;

  using matchingType = std::tuple<ttk::SimplexId, ttk::SimplexId, double>;

  double getPersistenceDiagram(std::vector<diagramType> &diagram,
                               vtkUnstructuredGrid *CTPersistenceDiagram_);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  void Modified() override;

  vtkSmartPointer<vtkUnstructuredGrid> createMatchings();
  vtkSmartPointer<vtkUnstructuredGrid> createOutputClusteredDiagrams();
  vtkSmartPointer<vtkUnstructuredGrid> createOutputCentroids();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::vector<std::vector<diagramType>> intermediateDiagrams_{};
  std::vector<std::vector<std::vector<matchingType>>> all_matchings_{};
  std::vector<std::vector<diagramType>> final_centroids_{};
  std::vector<int> inv_clustering_{};

  double Spacing{1.0};
  int DisplayMethod{0};
  double oldSpacing{1.0};

  double max_dimension_total_{};
  int Method{0}; // 0 = progressive approach, 1 = Auction approach
  bool needUpdate_{true};
};
