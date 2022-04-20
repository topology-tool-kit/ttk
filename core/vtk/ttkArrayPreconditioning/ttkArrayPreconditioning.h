/// \ingroup vtk
/// \class ttkArrayPreconditioning
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
/// \sa ttk::ttkArrayPreconditioning
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkArrayPreconditioningModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

#ifdef TTK_ENABLE_MPI
#include <mpi.h>
#endif
#include <vtkDataArraySelection.h>
#include <vtkNew.h>

class TTKARRAYPRECONDITIONING_EXPORT ttkArrayPreconditioning
  : public ttkAlgorithm {

public:
  static ttkArrayPreconditioning *New();
  vtkTypeMacro(ttkArrayPreconditioning, ttkAlgorithm);

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(BurstSize, int);
  vtkGetMacro(BurstSize, int);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

  // copy the vtkPassSelectedArray ("PassArrays" filter) API
  vtkDataArraySelection *GetPointDataArraySelection() {
    return this->ArraySelection;
  }

protected:
  ttkArrayPreconditioning();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

#ifdef TTK_ENABLE_MPI
  void ReceiveAndAddToVector(
    MPI_Datatype mpi_values,
    int rankFrom,
    int structTag,
    int intTag,
    std::vector<std::vector<ttk::value>> &unsortedReceivedValues);
#endif

private:
  vtkNew<vtkDataArraySelection> ArraySelection{};
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  int BurstSize{100};
};
