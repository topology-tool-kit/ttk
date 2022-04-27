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
/// The input data array that will be processed needs to be specified via the
/// standard VTK call SetInputArrayToProcess(), with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the VTK array)
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

  vtkGetMacro(RandomSeed, int);
  vtkSetMacro(RandomSeed, int);

  vtkGetMacro(CompactRange, bool);
  vtkSetMacro(CompactRange, bool);

protected:
  ttkIdentifierRandomizer();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int RandomSeed{};
  bool CompactRange{false};
};
