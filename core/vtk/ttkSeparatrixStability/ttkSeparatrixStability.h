/// \ingroup vtk
/// \class ttkSeparatrixStability
/// \author Thomas Daniel <thomas.daniel124@gmail.com>
/// \date July 2025.
///
/// \brief TTK VTK-filter that wraps the ttk::SeparatrixStability module.
///
/// It takes as an input a vtkMultiBlockDataSet representing the 1-dimensional 
/// separatrices of the Morse-Smale complex for an ensemble dataset. It 
/// computes as an output, a copy of the input, with for each separatrix, its 
/// rate of occurrence in the ensemble (based on partial isomorphism 
/// computations).
///
/// \sa ttk::SeparatrixStability
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/molecularVibration/">Molecular
///   Vibration example</a>
///   
/// \b Related \b publication: \n
/// "BondMatcher: H-Bond Stability Analysis in Molecular Systems" \n
/// Thomas Daniel, Malgorzata Olejniczak, Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics \n
/// Proc. of IEEE VIS 2025.

#pragma once

#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>

#include <ttkSeparatrixStabilityModule.h>

#include <SeparatrixStability.h>
#include <ttkAlgorithm.h>

class TTKSEPARATRIXSTABILITY_EXPORT ttkSeparatrixStability
  : public ttkAlgorithm,
    protected ttk::SeparatrixStability {

private:
  double PX{1};
  double PY{1};
  double PZ{1};
  double PF{1};
  bool MergeEdgesOnSaddles{true};
  double CostDeathBirth{};

public:
  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PF, double);
  vtkGetMacro(PF, double);

  vtkSetMacro(MergeEdgesOnSaddles, bool);
  vtkGetMacro(MergeEdgesOnSaddles, bool);

  vtkSetMacro(CostDeathBirth, double);
  vtkGetMacro(CostDeathBirth, double);

  static ttkSeparatrixStability *New();
  vtkTypeMacro(ttkSeparatrixStability, ttkAlgorithm);

protected:
  ttkSeparatrixStability();
  ~ttkSeparatrixStability() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int execute(vtkMultiBlockDataSet *&multiBlock1_Separatrices,
              vtkMultiBlockDataSet *&output1_Separatrices);

  bool updateVisitedVertices(const int &globalId,
                             std::vector<int> &localToGlobal,
                             int &localId);

  void computePointIds(const int &cellId_1,
                       const int &cellId_2,
                       const int &sourceGlobalId,
                       const int &destinationGlobalId,
                       vtkDataSet *block,
                       vtkIdType &srcPointId,
                       vtkIdType &destPointId);

  void updateAdjacencyMatrix(const int &sourceLocalId,
                             const int &destinationLocalId,
                             const int &separatrixLocalId,
                             GraphMatrixFull &adjacencyMatrix);

  void appendPoint(vtkPoints *points,
                   const int &index,
                   const double &scalar,
                   std::vector<std::array<double, 3>> &coords,
                   std::vector<double> &scalars);

  int prepareData(vtkDataSet *block,
                  std::vector<int> &localToGlobal,
                  GraphMatrixFull &adjacencyMatrixFull,
                  std::vector<std::array<double, 3>> &coordsSource,
                  std::vector<std::array<double, 3>> &coordsDestination,
                  std::vector<double> &scalarsSource,
                  std::vector<double> &scalarsDestinatoin,
                  int &n_separatrices,
                  std::vector<int> &globalSourcePointId,
                  std::vector<int> &globalDestinationPointId);
};
