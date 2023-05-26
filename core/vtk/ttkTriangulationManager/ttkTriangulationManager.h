/// \ingroup vtk
/// \class ttkTriangulationManager
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2022.
///
/// \brief TTK VTK-filter that manages ttk::Triangulation options.
///
///
/// \sa ttk::Triangulation
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/compactTriangulation/">
///   Compact Triangulation example</a> \n

#pragma once

// VTK Module
#include <ttkTriangulationManagerModule.h>

// VTK Includes
#include <vtkDataArraySelection.h>
#include <vtkSmartPointer.h>

#ifdef TTK_ENABLE_MPI
#include <ttkPeriodicGhostsGeneration.h>
#include <ttkTriangulationFactory.h>
#endif

// TTK Includes
#include <Triangulation.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class vtkPointSet;
class vtkUnstructuredGrid;

class TTKTRIANGULATIONMANAGER_EXPORT ttkTriangulationManager
  : public ttkAlgorithm {

  using STRATEGY = typename ttk::Triangulation::STRATEGY;
  bool Periodicity{false};
  STRATEGY PreconditioningStrategy{STRATEGY::DEFAULT};
  int Threshold{1000};
  vtkSmartPointer<vtkDataArraySelection> ArraySelection{};
#ifdef TTK_ENABLE_MPI
  vtkNew<ttkPeriodicGhostsGeneration> periodicGhostGenerator{};
#endif

public:
  ttkSetEnumMacro(PreconditioningStrategy, STRATEGY);
  vtkGetEnumMacro(PreconditioningStrategy, STRATEGY);

  vtkSetMacro(Periodicity, bool);
  vtkGetMacro(Periodicity, bool);

  vtkSetMacro(Threshold, int);
  vtkGetMacro(Threshold, int);

  // copy the vtkPassSelectedArray ("PassArrays" filter) API
  vtkDataArraySelection *GetDataArraySelection() {
    return this->ArraySelection.GetPointer();
  }

  void SetDataArraySelection(
    const vtkSmartPointer<vtkDataArraySelection> &selection) {
    this->ArraySelection = selection;
  }

  static ttkTriangulationManager *New();
  vtkTypeMacro(ttkTriangulationManager, ttkAlgorithm);

protected:
  ttkTriangulationManager();
  ~ttkTriangulationManager() override = default;

  void processImplicit(ttk::Triangulation &triangulation
#ifdef TTK_ENABLE_MPI
                       ,
                       vtkImageData *imageIn,
                       vtkImageData *imageOut
#endif
  );
  int processExplicit(vtkUnstructuredGrid *const output,
                      vtkPointSet *const input,
                      ttk::Triangulation &triangulation) const;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
#ifdef TTK_ENABLE_MPI
  int RequestUpdateExtent(vtkInformation *ttkNotUsed(request),
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector) override;
  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVectors,
                         vtkInformationVector *outputVector) override;
#endif // TTK_ENABLE_MPI
};
