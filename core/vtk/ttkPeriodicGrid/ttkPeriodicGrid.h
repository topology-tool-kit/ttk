///
/// \ingroup vtk
/// \class ttkPeriodicGrid
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date February 2020.
///
/// \brief TTK VTK-filter set the periodicity (in all dimensions) of a regular
/// grid.
///
///
/// \sa ttk::Triangulation
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkPeriodicGridModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkTriangulation.h>

class TTKPERIODICGRID_EXPORT ttkPeriodicGrid
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
{
private:
  bool Periodicity;

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */
  vtkSetMacro(Periodicity, bool);
  vtkGetMacro(Periodicity, bool);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkPeriodicGrid *New();
  vtkTypeMacro(ttkPeriodicGrid, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkPeriodicGrid();
  ~ttkPeriodicGrid() override;

  /**
   * TODO 8: Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 9: Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
