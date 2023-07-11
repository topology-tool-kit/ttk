/// \ingroup vtk
/// \class ttkIdentifiers
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date October 2022.
///
/// \brief TTK VTK-filter that triggers the computation of global identifiers.
///
/// This filter is useful when TTK is compiled with MPI to compute the
/// MPI preconditioning of the Identifiers filter.
/// When TTK is not compiled with MPI, the filter will compute the global
/// identifiers without triggering the MPI preconditioning.
///
/// \param Input Input data-set (vtkDataSet)
/// \param Output Output data-set with MPI preconditioning (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///

#pragma once

// VTK includes -- to adapt

// VTK Module
#include <ttkIdentifiersModule.h>

// ttk code includes
#include <Identifiers.h>
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.

class TTKIDENTIFIERS_EXPORT ttkIdentifiers : public ttkAlgorithm,
                                             protected ttk::Identifiers {

public:
  static ttkIdentifiers *New();

  vtkTypeMacro(ttkIdentifiers, ttkAlgorithm);

  vtkSetMacro(CellFieldName, const std::string &);
  vtkGetMacro(CellFieldName, std::string);

  vtkSetMacro(VertexFieldName, const std::string &);
  vtkGetMacro(VertexFieldName, std::string);

protected:
  ttkIdentifiers();

  ~ttkIdentifiers() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string CellFieldName{ttk::CellScalarFieldName},
    VertexFieldName{ttk::VertexScalarFieldName};
};
