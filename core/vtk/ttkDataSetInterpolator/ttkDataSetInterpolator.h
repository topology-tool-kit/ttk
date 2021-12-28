/// \ingroup vtk
/// \class ttkDataSetInterpolator
/// \author Guillaume Favelier <guillaume.favelier@gmail.com>
/// \date September 2018.
///
/// \brief TTK VTK-filter that wraps the dataSetInterpolator processing package.
///
/// VTK wrapping code for the ttk::DataSetInterpolator package.
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
/// \sa ttk::DataSetInterpolator
#pragma once

// VTK Module
#include <ttkDataSetInterpolatorModule.h>

// TTK Includes
#include <ttkAlgorithm.h>

class TTKDATASETINTERPOLATOR_EXPORT ttkDataSetInterpolator
  : public ttkAlgorithm {

public:
  static ttkDataSetInterpolator *New();
  vtkTypeMacro(ttkDataSetInterpolator, ttkAlgorithm);

protected:
  ttkDataSetInterpolator();
  ~ttkDataSetInterpolator();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
