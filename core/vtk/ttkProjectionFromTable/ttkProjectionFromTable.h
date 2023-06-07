/// \ingroup vtk
/// \class ttkProjectionFromTable
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// \brief TTK VTK-filter that wraps the
/// ttk::ProjectionFromTable module.
///
/// This VTK filter uses the ttk::ProjectionFromTable module
/// that projects on a surface points in a vtkTable.
///
/// \param Input vtkPolyData Surface
/// \param Input vtkTable Coefficients
/// \param Output vtkPolyData.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/ProjectionFromTable/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::ProjectionFromTable
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// VTK Module
#include <ttkProjectionFromTableModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkPolyData.h>
#include <vtkTable.h>

// TTK Base Includes
#include <ProjectionFromTable.h>

class TTKPROJECTIONFROMTABLE_EXPORT ttkProjectionFromTable
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ProjectionFromTable // and we inherit from
                                       // the base class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkProjectionFromTable *New();
  vtkTypeMacro(ttkProjectionFromTable, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkProjectionFromTable();
  ~ttkProjectionFromTable() override = default;

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
  template <typename VTK_T1, typename VTK_T2>
  void dispatch(std::vector<std::tuple<int, double, double>> &storage,
                const VTK_T1 *const values,
                const VTK_T2 *const values2,
                const size_t nvalues);
};
