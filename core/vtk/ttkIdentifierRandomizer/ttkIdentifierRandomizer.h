/// \ingroup vtk
/// \class ttkIdentifierRandomizer
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2017.
///
/// \brief TTK VTK-filter that randomly shuffles segmentation identifiers.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::IdentifierRandomizer
#pragma once

#include <ttkAlgorithm.h>
#include <ttkIdentifierRandomizerModule.h>

class TTKIDENTIFIERRANDOMIZER_EXPORT ttkIdentifierRandomizer
  : public ttkAlgorithm {

public:
  static ttkIdentifierRandomizer *New();
  vtkTypeMacro(ttkIdentifierRandomizer, ttkAlgorithm);

protected:
  ttkIdentifierRandomizer();
  ~ttkIdentifierRandomizer();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
