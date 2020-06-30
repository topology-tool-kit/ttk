/// \ingroup vtk
/// \class ttkIdentifyByScalarField
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that computes a new scalar array based on a sorting of the input array.
///
/// VTK code for computing a new scalar array from a specifiable point or cell data array. For a cell or a point
/// the new (integer) value is its index in a sorting of the input array.
/// The name of the new array is "PointScalarFieldName" or "CellScalarFieldName".
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// The name of the input scalar field to consider can be specified with the
/// standard VTK call SetInputArrayToProcess(), with the following parameters:
/// \param port 0
/// \param connection 0
/// \param fieldAssociation 0 or 1 (point data or cell data)
/// \param arrayName (const char* string representing the name of the VTK array)
///
/// This module respects the following convention regarding the order of the
/// input arrays to process (SetInputArrayToProcess()):
/// \param idx 0: input scalar field
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::IdentifyByScalarField
#pragma once

// VTK Module
#include <ttkIdentifyByScalarFieldModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKIDENTIFYBYSCALARFIELD_EXPORT ttkIdentifyByScalarField
  : public ttkAlgorithm {

public:
  static ttkIdentifyByScalarField *New();
  vtkTypeMacro(ttkIdentifyByScalarField, ttkAlgorithm)

  vtkSetMacro(IncreasingOrder, bool);
  vtkGetMacro(IncreasingOrder, bool);

  vtkSetMacro(StartByOne, bool);
  vtkGetMacro(StartByOne, bool);

  template <typename VTK_TT>
  int dispatch(std::vector<ttk::SimplexId> &inputIds);

protected:
  ttkIdentifyByScalarField();

  ~ttkIdentifyByScalarField() override{};

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,vtkInformationVector **inputVector,vtkInformationVector *outputVector) override;

private:
  bool IncreasingOrder{false};
  bool StartByOne{false};

  vtkDataArray *inputScalars_{nullptr};
};
