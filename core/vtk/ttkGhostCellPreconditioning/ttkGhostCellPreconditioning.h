/// \ingroup vtk
/// \class ttkGhostCellPreconditioning
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022
///
/// \brief TTK VTK-filter to generate rank ownership information.
///
/// This VTK filter generates a rankArray based on a global id array and ghost
/// cell information.
///
/// \param Input vtkDataSet with global ids and ghost cell information.
/// \param Output vtkDataSet.
///
///
/// \sa ttk::ttkGhostCellPreconditioning
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkGhostCellPreconditioningModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

#include <vtkDataArraySelection.h>
#include <vtkNew.h>

class TTKGHOSTCELLPRECONDITIONING_EXPORT ttkGhostCellPreconditioning
  : public ttkAlgorithm {

public:
  static ttkGhostCellPreconditioning *New();
  vtkTypeMacro(ttkGhostCellPreconditioning, ttkAlgorithm);

  // copy the vtkPassSelectedArray ("PassArrays" filter) API
  vtkDataArraySelection *GetPointDataArraySelection() {
    return this->ArraySelection;
  }

protected:
  ttkGhostCellPreconditioning();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  vtkNew<vtkDataArraySelection> ArraySelection{};
};
