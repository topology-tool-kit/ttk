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
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceDiagramDistanceMatrixModule.h>

// ttk code includes
#include <PersistenceDiagramDistanceMatrix.h>
#include <ttkAlgorithm.h>

class TTKPERSISTENCEDIAGRAMDISTANCEMATRIX_EXPORT
  ttkPersistenceDiagramDistanceMatrix
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramDistanceMatrix {

public:
  static ttkPersistenceDiagramDistanceMatrix *New();

  vtkTypeMacro(ttkPersistenceDiagramDistanceMatrix, ttkAlgorithm);

  vtkSetMacro(WassersteinMetric, std::string);
  vtkGetMacro(WassersteinMetric, std::string);

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

  vtkSetMacro(PairTypeClustering, int);
  vtkGetMacro(PairTypeClustering, int);

  vtkSetMacro(UseFullDiagrams, bool);
  vtkGetMacro(UseFullDiagrams, bool);

protected:
  ttkPersistenceDiagramDistanceMatrix();
  ~ttkPersistenceDiagramDistanceMatrix() override = default;

  double getPersistenceDiagram(std::vector<ttk::DiagramTuple> &diagram,
                               vtkUnstructuredGrid *CTPersistenceDiagram_);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int PairTypeClustering{-1};
  double Alpha{1.0};
  double DeltaLim{0.01};
  double Lambda{1.0};
  std::string WassersteinMetric{"2"};
  bool UseFullDiagrams{false};
};
