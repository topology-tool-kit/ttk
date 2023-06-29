/// \ingroup vtk
/// \class ttkPersistenceDiagramDictionary
/// \author Keanu Sisouk <keanu.sisouk@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date Mai 2023
///
/// \brief TTK processing package for the computation of a Dictionary
/// of Persistence Diagrams and barycentric weights to approximate
/// an ensemble of Persistence Diagrams.
///
/// \b Related \b publication \n
/// "Wasserstein Dictionaries of Persistence Diagrams" \n
/// Keanu Sisouk, Julie Delon and Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa PersistenceDiagramDictionary

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkSetGet.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceDiagramDictionaryModule.h>

// ttk code includes

#include <PersistenceDiagramDictionary.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class TTKPERSISTENCEDIAGRAMDICTIONARY_EXPORT ttkPersistenceDiagramDictionary
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramDictionary {

private:
  int AtomNumber_{3};
  int Seed_{0};
  double Percent_{0};

public:
  // enum class BACKEND{BORDER_INIT = 0 , RANDOM_INIT = 1 , FIRST_DIAGS = 2};

  static ttkPersistenceDiagramDictionary *New();

  vtkTypeMacro(ttkPersistenceDiagramDictionary, ttkAlgorithm);

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

  vtkSetMacro(Percent_, double);
  vtkGetMacro(Percent_, double);

  vtkSetMacro(OptimizeWeights_, int);
  vtkGetMacro(OptimizeWeights_, int);

  vtkSetMacro(OptimizeAtoms_, int);
  vtkGetMacro(OptimizeAtoms_, int);

  vtkSetMacro(MaxEigenValue_, int);
  vtkGetMacro(MaxEigenValue_, int);

  vtkSetMacro(Fusion_, int);
  vtkGetMacro(Fusion_, int);

  vtkSetMacro(ProgBarycenter_, int);
  vtkGetMacro(ProgBarycenter_, int);

  vtkSetMacro(MaxEpoch_, int);
  vtkGetMacro(MaxEpoch_, int);

  vtkSetMacro(ProgApproach_, int);
  vtkGetMacro(ProgApproach_, int);

  vtkSetMacro(StopCondition_, int);
  vtkGetMacro(StopCondition_, int);

  vtkSetMacro(CompressionMode_, int);
  vtkGetMacro(CompressionMode_, int);

  vtkSetMacro(DimReductMode_, int);
  vtkGetMacro(DimReductMode_, int);

  vtkSetMacro(sortedForTest_, int);
  vtkGetMacro(sortedForTest_, int);

  vtkSetMacro(CreationFeatures_, int);
  vtkGetMacro(CreationFeatures_, int);

  vtkSetMacro(AtomNumber_, int);
  vtkGetMacro(AtomNumber_, int);

  vtkSetMacro(Seed_, int);
  vtkGetMacro(Seed_, int);

  vtkSetMacro(DeltaLim, double);
  vtkGetMacro(DeltaLim, double);

  vtkSetMacro(Lambda, double);
  vtkGetMacro(Lambda, double);

  ttkSetEnumMacro(BackEnd, BACKEND);
  vtkGetEnumMacro(BackEnd, BACKEND);

  vtkSetMacro(CompressionFactor, double);
  vtkGetMacro(CompressionFactor, double);

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

  vtkSetMacro(MaxNumberOfPairs, unsigned int);
  vtkGetMacro(MaxNumberOfPairs, unsigned int);

  vtkSetMacro(MinPersistence_, double);
  vtkGetMacro(MinPersistence_, double);

protected:
  ttkPersistenceDiagramDictionary();
  ~ttkPersistenceDiagramDictionary() override = default;

  // BACKEND BackEnd{BACKEND::BORDER_INIT};
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
