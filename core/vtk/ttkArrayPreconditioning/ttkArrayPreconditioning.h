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

#include <vtkDataArraySelection.h>
#include <vtkNew.h>

// TTK Base Includes
#include <ArrayPreconditioning.h>

class TTKARRAYPRECONDITIONING_EXPORT ttkArrayPreconditioning
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ArrayPreconditioning // and we inherit from the base class
{

public:
  static ttkArrayPreconditioning *New();
  vtkTypeMacro(ttkArrayPreconditioning, ttkAlgorithm);

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(BurstSize, int);
  vtkGetMacro(BurstSize, int);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

  vtkSetMacro(GlobalOrder, bool);
  vtkGetMacro(GlobalOrder, bool);

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

private:
  vtkNew<vtkDataArraySelection> ArraySelection{};
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  int BurstSize{100};
};
