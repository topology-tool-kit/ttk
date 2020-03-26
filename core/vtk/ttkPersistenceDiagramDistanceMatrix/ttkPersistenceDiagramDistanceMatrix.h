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
    needUpdate_ = true;
  }

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
    needUpdate_ = true;
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

  void SetTimeLimit(double data) {
    TimeLimit = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(TimeLimit, double);

  void SetThreadNumber(int data) {
    ThreadNumber = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(ThreadNumber, int);

  void SetAlpha(double data) {
    if(data > 0 && data <= 1) {

      Alpha = data;
    } else if(data > 1) {
      Alpha = 1;
    } else {
      Alpha = 0.001;
    }
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Alpha, double);

  void SetAntiAlpha(double data) {
    double alpha = 1 - data;
    SetAlpha(alpha);
  }

  void SetDeltaLim(double data) {
    DeltaLim = data;
    Modified();
    needUpdate_ = true;
  }

  vtkGetMacro(DeltaLim, double);

  void SetLambda(double data) {
    Lambda = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Lambda, double);

  void SetUseAccelerated(bool data) {
    UseAccelerated = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseAccelerated, bool);

  void SetUseKmeansppInit(bool data) {
    UseKmeansppInit = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseKmeansppInit, bool);

  void SetForceUseOfAlgorithm(bool data) {
    ForceUseOfAlgorithm = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(ForceUseOfAlgorithm, bool);

  void SetDeterministic(bool data) {
    Deterministic = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Deterministic, bool);

  void SetPairTypeClustering(int data) {
    PairTypeClustering = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(PairTypeClustering, int);

  void SetUseAdditionalPrecision(bool data) {
    UseAdditionalPrecision = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseAdditionalPrecision, bool);

  void SetDistanceWritingOptions(int data) {
    DistanceWritingOptions = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(DistanceWritingOptions, int);

  void SetUseInterruptible(bool data) {
    UseInterruptible = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseInterruptible, bool);

  void SetUseFullDiagrams(const bool arg) {
    UseFullDiagrams = arg;
    needUpdate_ = true;
    Modified();
  }
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
  std::vector<std::vector<ttk::DiagramTuple>> intermediateDiagrams_{};

  int numberOfInputsFromCommandLine{1};
  int PairTypeClustering{-1};
  bool ForceUseOfAlgorithm{false};
  bool Deterministic{true};
  bool UseAllCores{false};
  int ThreadNumber{1};
  bool UseAdditionalPrecision{false};
  int DistanceWritingOptions{0};
  double Alpha{1.0};
  double DeltaLim{0.01};
  double Lambda{1.0};
  bool UseInterruptible{true};
  double max_dimension_total_{};

  bool needUpdate_{true};

  bool UseAccelerated{false};
  bool UseKmeansppInit{false};
  std::string WassersteinMetric{"2"};

  double TimeLimit{9999999};
  bool UseFullDiagrams{false};
  std::vector<std::vector<double>> diagramsDistMat{};

  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};
