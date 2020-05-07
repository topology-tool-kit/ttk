/// \ingroup vtk
/// \class ttkIdentifyByScalarField
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the identifyByScalarField processing
/// package.
///
/// VTK wrapping code for the @IdentifyByScalarField package.
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
/// \sa ttk::IdentifyByScalarField
#pragma once

#include <ttkAlgorithm.h>
#include <ttkIdentifyByScalarFieldModule.h>

class TTKIDENTIFYBYSCALARFIELD_EXPORT ttkIdentifyByScalarField
  : public ttkAlgorithm {

private:
  int ScalarFieldId{0};
  bool IncreasingOrder{false};
  bool StartByOne{false};

  vtkDataArray *inputScalars_;

public:
  vtkSetMacro(IncreasingOrder, bool);
  vtkGetMacro(IncreasingOrder, bool);

  vtkSetMacro(StartByOne, bool);
  vtkGetMacro(StartByOne, bool);

  static ttkIdentifyByScalarField *New();
  vtkTypeMacro(ttkIdentifyByScalarField, ttkAlgorithm);

protected:
  ttkIdentifyByScalarField();
  ~ttkIdentifyByScalarField();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
