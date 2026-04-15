/// \ingroup vtk
/// \class ttkPersistenceDiagramDictionaryDecoding
/// \author Keanu Sisouk <keanu.sisouk@lip6.fr>
/// \date February 2026
///
/// \brief TTK processing package for the computation of a Dictionary
/// of Persistence Diagrams and barycentric weights to approximate
/// an ensemble of Persistence Diagrams.
///
/// \b Related \b publication \n
/// "Wasserstein Dictionaries of Persistence Diagrams" \n
/// Keanu Sisouk, Julie Delon and Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics, 2024.
///
/// \sa PersistenceDiagramDictionaryDecoding
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramDictionary/">Persistence
///   Diagram Dictionary example</a> \n

#pragma once

#include "ttkMacros.h"
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>
// VTK Module
#include <ttkPersistenceDiagramDictionaryDecodingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <PersistenceDiagramDictionaryDecoding.h>

class TTKPERSISTENCEDIAGRAMDICTIONARYDECODING_EXPORT
  ttkPersistenceDiagramDictionaryDecoding
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramDictionaryDecoding {

private:
public:
  static ttkPersistenceDiagramDictionaryDecoding *New();
  vtkTypeMacro(ttkPersistenceDiagramDictionaryDecoding, ttkAlgorithm);

  vtkGetMacro(Spacing, double);
  vtkSetMacro(Spacing, double);

  vtkGetMacro(ShowAtoms, int);
  vtkSetMacro(ShowAtoms, int);

  vtkGetMacro(ProgBarycenter, int);
  vtkSetMacro(ProgBarycenter, int);

  ttkSetEnumMacro(ProjMet, BACKEND);
  vtkGetEnumMacro(ProjMet, BACKEND);

protected:
  ttkPersistenceDiagramDictionaryDecoding();
  ~ttkPersistenceDiagramDictionaryDecoding() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  void GetOutputDiagrams(vtkMultiBlockDataSet *output,
                         vtkTable *output_coordinates,
                         const std::vector<ttk::DiagramType> &diags,
                         std::vector<ttk::DiagramType> &atoms,
                         vtkTable *weights_vtk,
                         const std::vector<std::vector<double>> &weights,
                         const double spacing,
                         const double maxPersistence) const;

  // double GetPersistenceOfGlobalPair(const ttk::DiagramType &diagram) const;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  double Spacing{};
  int ShowAtoms{1};
  bool ComputePoints{false};
};
