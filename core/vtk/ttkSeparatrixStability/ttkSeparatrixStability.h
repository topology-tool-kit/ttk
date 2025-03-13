/// \ingroup vtk
/// \class ttkSeparatrixStability
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::SeparatrixStability module.
///
/// This VTK filter uses the ttk::SeparatrixStability module to compute
/// an averaging of the data values of an input point data array defined on the
/// input vtkDataSet.
///
/// \sa ttk::SeparatrixStability
/// \sa ttkAlgorithm

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
  bool ComputeOccurenceType0{true};
  bool ComputeOccurenceType1{false};
  bool ComputeOccurenceType2{true};
  bool MergeEdgesOnSaddles{true};
  double CostDeathBirth{};
  bool ShuffleCriticalPointsIds{true};
  bool ShuffleSeparatrixIds{true};
  unsigned int SeparatrixShuffleSeed{0};
  unsigned int CriticalPointsShuffleSeed{0};

public:

  vtkSetMacro(MergeEdgesOnSaddles, bool);
  vtkGetMacro(MergeEdgesOnSaddles, bool);

  vtkSetMacro(CostDeathBirth, double);
  vtkGetMacro(CostDeathBirth, double);

  vtkSetMacro(ShuffleCriticalPointsIds, bool);
  vtkGetMacro(ShuffleCriticalPointsIds, bool);

  vtkSetMacro(CriticalPointsShuffleSeed, unsigned int);
  vtkGetMacro(CriticalPointsShuffleSeed, unsigned int);

  vtkSetMacro(ShuffleSeparatrixIds, bool);
  vtkGetMacro(ShuffleSeparatrixIds, bool);

  vtkSetMacro(SeparatrixShuffleSeed, unsigned int);
  vtkGetMacro(SeparatrixShuffleSeed, unsigned int);

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
                   std::vector<std::array<double, 3>> &coords);

  int prepareData(vtkDataSet *block,
                  std::vector<int> &localToGlobal,
                  GraphMatrixFull &adjacencyMatrixFull,
                  std::vector<std::array<double, 3>> &coordsSource,
                  std::vector<std::array<double, 3>> &coordsDestination,
                  int &n_separatrices,
                  std::vector<int> &globalSourcePointId,
                  std::vector<int> &globalDestinationPointId);
};
