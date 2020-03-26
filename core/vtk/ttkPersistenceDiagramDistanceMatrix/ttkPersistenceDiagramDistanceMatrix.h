/// \ingroup base
/// \class ttkPersistenceDiagramDistanceMatrix
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
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
/// \sa PersistenceDiagramDistanceMatrix

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceDiagramDistanceMatrixModule.h>

// ttk code includes
#include <PersistenceDiagramDistanceMatrix.h>
#include <ttkTriangulationAlgorithm.h>

class TTKPERSISTENCEDIAGRAMDISTANCEMATRIX_EXPORT
  ttkPersistenceDiagramDistanceMatrix : public vtkMultiBlockDataSetAlgorithm,
                                        protected ttk::Wrapper {

public:
  static ttkPersistenceDiagramDistanceMatrix *New();

  vtkTypeMacro(ttkPersistenceDiagramDistanceMatrix,
               vtkMultiBlockDataSetAlgorithm);

  // default ttk setters
  void SetDebugLevel(int debugLevel) {
    setDebugLevel(debugLevel);
    Modified();
  }

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
  }

  /*void SetThreadNumber(int threadNumber){
    ThreadNumber = threadNumber;
    SetThreads();
  }*/

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  // set-getters macros to define from each variable you want to access from
  // the outside (in particular from paraview) - to adapt.

  vtkSetMacro(WassersteinMetric, std::string);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(TimeLimit, double);
  vtkGetMacro(TimeLimit, double);

  void SetThreadNumber(int data) {
    ThreadNumber = data;
    Modified();
  }
  vtkGetMacro(ThreadNumber, int);

  void SetAntiAlpha(double data) {
    data = 1 - data;
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

  vtkSetMacro(DeltaLim, double);
  vtkGetMacro(DeltaLim, double);

  vtkSetMacro(Lambda, double);
  vtkGetMacro(Lambda, double);

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

  vtkSetMacro(UseAdditionalPrecision, bool);
  vtkGetMacro(UseAdditionalPrecision, bool);

  vtkSetMacro(UseInterruptible, bool);
  vtkGetMacro(UseInterruptible, bool);

  vtkSetMacro(UseFullDiagrams, bool);
  vtkGetMacro(UseFullDiagrams, bool);

protected:
  ttkPersistenceDiagramDistanceMatrix();

  double getPersistenceDiagram(std::vector<ttk::DiagramTuple> &diagram,
                               vtkUnstructuredGrid *CTPersistenceDiagram_);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int PairTypeClustering{-1};
  bool ForceUseOfAlgorithm{false};
  bool Deterministic{true};
  bool UseAllCores{false};
  int ThreadNumber{1};
  bool UseAdditionalPrecision{false};
  double Alpha{1.0};
  double DeltaLim{0.01};
  double Lambda{1.0};
  bool UseInterruptible{true};

  bool UseAccelerated{false};
  bool UseKmeansppInit{false};
  std::string WassersteinMetric{"2"};

  double TimeLimit{9999999};
  bool UseFullDiagrams{false};

  bool needsToAbort() override;
  int updateProgress(const float &progress) override;
};
