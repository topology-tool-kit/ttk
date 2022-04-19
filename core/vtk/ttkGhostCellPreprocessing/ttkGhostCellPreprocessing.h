/// \ingroup vtk
/// \class ttkGhostCellPreprocessing
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date September 2020
///
/// \brief TTK VTK-filter to generate order fields.
///
/// This VTK filter generates order arrays from a selection of scalar
/// field arrays.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::ttkGhostCellPreprocessing
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkGhostCellPreprocessingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

#include <mpi.h>
#include <vtkDataArraySelection.h>
#include <vtkNew.h>

class TTKGHOSTCELLPREPROCESSING_EXPORT ttkGhostCellPreprocessing
  : public ttkAlgorithm {

public:
  static ttkGhostCellPreprocessing *New();
  vtkTypeMacro(ttkGhostCellPreprocessing, ttkAlgorithm);

  // copy the vtkPassSelectedArray ("PassArrays" filter) API
  vtkDataArraySelection *GetPointDataArraySelection() {
    return this->ArraySelection;
  }

protected:
  ttkGhostCellPreprocessing();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  vtkNew<vtkDataArraySelection> ArraySelection{};
};
