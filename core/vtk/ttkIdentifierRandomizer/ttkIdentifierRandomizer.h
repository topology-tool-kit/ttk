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

// VTK includes -- to adapt

// VTK Module
#include <ttkIdentifierRandomizerModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKIDENTIFIERRANDOMIZER_EXPORT ttkIdentifierRandomizer
  : public ttkAlgorithm {

public:
  static ttkIdentifierRandomizer *New();
  vtkTypeMacro(ttkIdentifierRandomizer, ttkAlgorithm);

protected:
  ttkIdentifierRandomizer();

  ~ttkIdentifierRandomizer() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
};
