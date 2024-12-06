/// \ingroup vtk
/// \class ttkMorseSmallComplexStability
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::MorseSmallComplexStability module.
///
/// This VTK filter uses the ttk::MorseSmallComplexStability module to compute an averaging of
/// the data values of an input point data array defined on the input
/// vtkDataSet.
///
/// \sa ttk::MorseSmallComplexStability
/// \sa ttkAlgorithm

#pragma once

#include<vtkMultiBlockDataSet.h>
#include<vtkUnstructuredGrid.h>

#include <ttkMorseSmallComplexStabilityModule.h>

#include <ttkAlgorithm.h>
#include <MorseSmallComplexStability.h>

class TTKMORSESMALLCOMPLEXSTABILITY_EXPORT ttkMorseSmallComplexStability
  : public ttkAlgorithm, 
    protected ttk::MorseSmallComplexStability // and we inherit from the base class
{
private:
  bool ComputeOccurenceType0{true};
  bool ComputeOccurenceType1{false};
  bool ComputeOccurenceType2{true};
  bool MergeEdgesOnSaddles{true};

public:
  //vtkSetMacro(ComputeOccurenceType0, bool);
  //vtkGetMacro(ComputeOccurenceType0, bool);
//
  //vtkSetMacro(ComputeOccurenceType1, bool);
  //vtkGetMacro(ComputeOccurenceType1, bool);
//
  //vtkSetMacro(ComputeOccurenceType2, bool);
  //vtkGetMacro(ComputeOccurenceType2, bool);

  vtkSetMacro(MergeEdgesOnSaddles, bool);
  vtkGetMacro(MergeEdgesOnSaddles, bool);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMorseSmallComplexStability *New();
  vtkTypeMacro(ttkMorseSmallComplexStability, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkMorseSmallComplexStability();
  ~ttkMorseSmallComplexStability() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;


  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int execute(vtkMultiBlockDataSet* &multiBlock1_Separatrices,
              vtkMultiBlockDataSet* &output1_Separatrices);
  
  bool updateVisitedVertices(const int &globalId, 
                              std::vector<int> &localToGlobal,
                              int &localId);

  void updateVertexData(const int &globalId, 
                        const vtkidType &pointId, 
                        const vtkPoints* points,
                        std::vector<std::array<double, 3>> &coords, 
                        std::vector<int> &localToGlobal,
                        int &localId);

  void computePointIds(const vtkCell* cell_1,
                      const vtkCell* cell_2,
                      const int &sourceGlobalId,
                      const int &destinationGlobalId,
                      const vtkDataSet &block,
                      vtkIdType &srcPointId,
                      vtkIdType &destPointId);

  void updateAdjacencyMatrix(const int &sourceLocalId,
                              const int &destinationLocalId,
                              const int &separatrixLocalId,
                              GraphMatrixFull &adjacencyMatrix);
                                                              
  void appendPoint(vtkPoints* points, 
                    const int &index, 
                    std::vector<std::array<double, 3>> &coords);

  void computeGraphMinor(const GraphMatrixFull &adjacencyMatrixFull, 
                          GraphMatrixMinor &adjacencyMatrix);

  int prepareData(vtkDataSet* block, 
                  std::vector<int> &localToGlobal, 
                  GraphMatrixFull &adjacencyMatrixFull,
                  GraphMatrixMinor &adjacencyMatrixMinor,
                  std::vector<std::array<double, 3>> &coordsSource,
                  std::vector<std::array<double, 3>> &coordsDestination,
                  int &n_separatrices);
  
};
