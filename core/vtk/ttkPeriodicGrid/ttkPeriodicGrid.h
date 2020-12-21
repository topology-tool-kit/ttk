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

class TTKPERIODICGRID_EXPORT ttkPeriodicGrid
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
{
private:
  bool Periodicity{true};

public:
  vtkSetMacro(Periodicity, bool);
  vtkGetMacro(Periodicity, bool);

  static ttkPeriodicGrid *New();
  vtkTypeMacro(ttkPeriodicGrid, ttkAlgorithm);

protected:
  ttkPeriodicGrid();
  ~ttkPeriodicGrid() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
