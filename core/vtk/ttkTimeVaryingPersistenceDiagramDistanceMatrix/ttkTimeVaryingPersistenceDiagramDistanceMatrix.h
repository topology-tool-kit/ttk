/// \ingroup vtk
/// \class ttkTimeVaryingPersistenceDiagramDistanceMatrix
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief TTK processing package for the computation of a matrix of Wasserstein
/// distances between persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa TimeVaryingPersistenceDiagramDistanceMatrix
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkInformationVector.h>

// VTK Module
#include <ttkTimeVaryingPersistenceDiagramDistanceMatrixModule.h>

// ttk code includes
#include <TimeVaryingPersistenceDiagramDistanceMatrix.h>
#include <ttkAlgorithm.h>

class TTKTIMEVARYINGPERSISTENCEDIAGRAMDISTANCEMATRIX_EXPORT
  ttkTimeVaryingPersistenceDiagramDistanceMatrix
  : public ttkAlgorithm,
    protected ttk::TimeVaryingPersistenceDiagramDistanceMatrix {

public:
  static ttkTimeVaryingPersistenceDiagramDistanceMatrix *New();

  vtkTypeMacro(ttkTimeVaryingPersistenceDiagramDistanceMatrix, ttkAlgorithm);

  void SetWassersteinMetric(const std::string &data) {
    Wasserstein = (data == "inf") ? -1 : stoi(data);
    Modified();
  }
  std::string GetWassersteinMetric() {
    return Wasserstein == -1 ? "inf" : std::to_string(Wasserstein);
  }

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

  void SetPairType(const int data) {
    switch(data) {
      case(0):
        this->setDos(true, false, false);
        break;
      case(1):
        this->setDos(false, true, false);
        break;
      case(2):
        this->setDos(false, false, true);
        break;
      default:
        this->setDos(true, true, true);
        break;
    }
    Modified();
  }
  int GetPairType() {
    if(do_min_ && do_sad_ && do_max_) {
      return -1;
    } else if(do_min_) {
      return 0;
    } else if(do_sad_) {
      return 1;
    } else if(do_max_) {
      return 2;
    }
    return -1;
  }

  void SetConstraint(const int arg_) {
    this->setConstraint(arg_);
    this->Modified();
  }
  int GetConstraint() {
    switch(this->Constraint) {
      case ConstraintType::FULL_DIAGRAMS:
        return 0;
      case ConstraintType::NUMBER_PAIRS:
        return 1;
      case ConstraintType::ABSOLUTE_PERSISTENCE:
        return 2;
      case ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG:
        return 3;
      case ConstraintType::RELATIVE_PERSISTENCE_GLOBAL:
        return 4;
    }
    return -1;
  }

  vtkSetMacro(MaxNumberOfPairs, unsigned int);
  vtkGetMacro(MaxNumberOfPairs, unsigned int);

  vtkSetMacro(MinPersistence, double);
  vtkGetMacro(MinPersistence, double);

protected:
  ttkTimeVaryingPersistenceDiagramDistanceMatrix();
  ~ttkTimeVaryingPersistenceDiagramDistanceMatrix() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
