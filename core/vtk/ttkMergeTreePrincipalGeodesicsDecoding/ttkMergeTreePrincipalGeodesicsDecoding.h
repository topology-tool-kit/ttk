/// \ingroup vtk
/// \class ttkMergeTreePrincipalGeodesicsDecoding
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// \brief TTK VTK-filter that wraps the
/// ttk::MergeTreePrincipalGeodesicsDecoding module.
///
/// This VTK filter uses the ttk::MergeTreePrincipalGeodesicsDecoding module to
/// compute the reconstruction of merge trees or persistence diagrams given a
/// Principal Geodesic Analysis basis and projection coefficients.
///
/// \param Input vtkMultiBlockDataSet Barycenter
/// \param Input vtkMultiBlockDataSet Geodesics
/// \param Input vtkTable Coefficients
/// \param Input vtkTable Branches Correlation
/// \param Output vtkMultiBlockDataSet Trees
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MergeTreePrincipalGeodesicsDecoding/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreePrincipalGeodesicsDecoding
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreePrincipalGeodesicsDecodingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <MergeTreePrincipalGeodesicsDecoding.h>
#include <ttkMergeTreeVisualization.h>

class TTKMERGETREEPRINCIPALGEODESICSDECODING_EXPORT
  ttkMergeTreePrincipalGeodesicsDecoding
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreePrincipalGeodesicsDecoding // and we inherit from
                                                       // the base class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */

  // Output options
  bool OutputInputTrees = false;
  bool OutputInputTreesSegmentation = false;
  bool OutputBarycenter = true;
  bool ReconstructInputTrees = true;
  bool ConstructGeodesicsTrees = false;
  bool ConstructEllipses = false;
  bool ConstructRectangle = false;
  unsigned int RectangleMultiplier = 1;
  bool ConstructSurface = false;
  bool ProcessSecondInput = false;

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<ttk::ftm::MergeTree<double>> baryMTree, inputMTrees;
  std::vector<vtkUnstructuredGrid *> baryTreeNodes, inputTreesNodes;
  std::vector<vtkUnstructuredGrid *> baryTreeArcs, inputTreesArcs;
  std::vector<vtkDataSet *> baryTreeSegmentation, inputTreesSegmentation;
  // Output
  bool processFirstInput;
  std::vector<ttk::ftm::MergeTree<double>> reconstructedTrees;
  std::vector<double> reconstructionErrors;
  std::vector<
    std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>
    recInputMatchings, recBaryMatchings;
  std::vector<std::vector<ttk::ftm::MergeTree<double>>> geodesicsTrees;
  std::vector<ttk::ftm::MergeTree<double>> geodesicsEllipses,
    geodesicsRectangle, geodesicsSurface;
  vtkFieldData *inputFieldData;
  // Verify input changed
  vtkMTimeType tableCoefficientsMTime, tableVectorsMTime, tableCorrelationMTime,
    blockInputTreesMTime;

  void setDataVisualization(int numTrees) {
    // Trees
    baryMTree.resize(1);
    baryTreeNodes.resize(1);
    baryTreeArcs.resize(1);
    baryTreeSegmentation.resize(1);

    inputMTrees.resize(numTrees);
    inputTreesNodes.resize(numTrees);
    inputTreesArcs.resize(numTrees);
    inputTreesSegmentation.resize(numTrees);
  }

  bool isDataVisualizationFilled() {
    return (ReconstructInputTrees == (reconstructedTrees.size() != 0))
           and (ConstructGeodesicsTrees == (geodesicsTrees.size() != 0))
           and (ConstructEllipses == (geodesicsEllipses.size() != 0))
           and (ConstructRectangle == (geodesicsRectangle.size() != 0))
           and (ConstructSurface == (geodesicsSurface.size() != 0))
           and (not OutputInputTrees) and geodesicsDistances_.size() != 0;
  }

  void resetDataVisualization() {
    setDataVisualization(0);
    reconstructedTrees.clear();
    geodesicsTrees.clear();
    geodesicsEllipses.clear();
    geodesicsRectangle.clear();
    geodesicsSurface.clear();
    geodesicsDistances_.clear();
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input options
  void SetNumberOfGeodesicsIntervals(int k) {
    k_ = k;
    Modified();
    resetDataVisualization();
  }

  int GetNumberOfGeodesicsIntervals() {
    return k_;
  }

  // Output options
  vtkSetMacro(OutputInputTrees, bool);
  vtkGetMacro(OutputInputTrees, bool);

  vtkSetMacro(OutputInputTreesSegmentation, bool);
  vtkGetMacro(OutputInputTreesSegmentation, bool);

  vtkSetMacro(OutputBarycenter, bool);
  vtkGetMacro(OutputBarycenter, bool);

  vtkSetMacro(ReconstructInputTrees, bool);
  vtkGetMacro(ReconstructInputTrees, bool);

  void SetcomputeReconstructionError_(bool b) {
    computeReconstructionError_ = b;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(computeReconstructionError_, bool);

  void SettransferInputTreesInformation_(bool b) {
    transferInputTreesInformation_ = b;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(transferInputTreesInformation_, bool);

  void SettransferBarycenterInformation_(bool b) {
    transferBarycenterInformation_ = b;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(transferBarycenterInformation_, bool);

  vtkSetMacro(ConstructGeodesicsTrees, bool);
  vtkGetMacro(ConstructGeodesicsTrees, bool);

  vtkSetMacro(ConstructEllipses, bool);
  vtkGetMacro(ConstructEllipses, bool);

  vtkSetMacro(ConstructRectangle, bool);
  vtkGetMacro(ConstructRectangle, bool);

  void SetRectangleMultiplier(unsigned int mult) {
    RectangleMultiplier = mult;
    Modified();
    geodesicsRectangle.clear();
  }
  vtkGetMacro(RectangleMultiplier, unsigned int);

  vtkSetMacro(ConstructSurface, bool);
  vtkGetMacro(ConstructSurface, bool);

  void SetProcessSecondInput(bool b) {
    ProcessSecondInput = b;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(ProcessSecondInput, bool);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreePrincipalGeodesicsDecoding *New();
  vtkTypeMacro(ttkMergeTreePrincipalGeodesicsDecoding, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreePrincipalGeodesicsDecoding();

  /**
   * Specify the input data type of each input port
   * (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   * (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   * (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class dataType>
  int run(vtkInformationVector *outputVector,
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputBary,
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees);

  template <class dataType>
  int runCompute(
    vtkInformationVector *outputVector,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputBary,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees);

  template <class dataType>
  int runOutput(vtkInformationVector *outputVector,
                std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputBary,
                std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees);
};
