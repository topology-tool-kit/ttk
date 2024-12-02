/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
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

// VTK Module
#include <ttkMorseSmallComplexStabilityModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

/* Note on including VTK modules
 *
 * Each VTK module that you include a header from needs to be specified in this
 * module's vtk.module file, either in the DEPENDS or PRIVATE_DEPENDS (if the
 * header is included in the cpp file only) sections.
 *
 * In order to find the corresponding module, check its location within the VTK
 * source code. The VTK module name is composed of the path to the header. You
 * can also find the module name within the vtk.module file located in the same
 * directory as the header file.
 *
 * For example, vtkSphereSource.h is located in directory VTK/Filters/Sources/,
 * so its corresponding VTK module is called VTK::FiltersSources. In this case,
 * the vtk.module file would need to be extended to
 *
 * NAME
 *   ttkMorseSmallComplexStability
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <MorseSmallComplexStability.h>

class TTKMORSESMALLCOMPLEXSTABILITY_EXPORT ttkMorseSmallComplexStability
  : public ttkAlgorithm, 
    protected ttk::MorseSmallComplexStability // and we inherit from the base class
{
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

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

  using GraphMatrix = std::vector<std::vector<std::optional<std::pair<int, int>>>>;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int execute(vtkMultiBlockDataSet* &multiBlock1_Separatrices);
  
  void updateVisitedVertices(const int &globalId, 
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

  void updateVertexLinksList(std::vector<std::vector<int>> &vertexLinks, 
                              const int &v1, 
                              const int &v2);

  int prepareData(vtkDataSet* block, 
                  std::vector<int> &localToGlobal, 
                  std::vector<std::pair<int, int>> &edges,
                  std::vector<std::array<double, 3>> &coords,
                  std::vector<float>&sfValues);
};
