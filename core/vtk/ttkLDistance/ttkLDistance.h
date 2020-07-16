/// \ingroup vtk
/// \class ttkLDistance
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 26/02/2017
///
/// \brief TTK VTK-filter that wraps the lDistance processing package.
///
/// VTK wrapping code for the @LDistance package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::LDistance
#pragma once

// VTK Module
#include <ttkLDistanceModule.h>

// ttk code includes
#include <LDistance.h>
#include <ttkAlgorithm.h>

class TTKLDISTANCE_EXPORT ttkLDistance : public ttkAlgorithm,
                                         protected ttk::LDistance {
public:
  static ttkLDistance *New();
  vtkTypeMacro(ttkLDistance, ttkAlgorithm);

  vtkSetMacro(ScalarField1, std::string);
  vtkGetMacro(ScalarField1, std::string);

  vtkSetMacro(ScalarField2, std::string);
  vtkGetMacro(ScalarField2, std::string);

  vtkSetMacro(ScalarFieldId1, int);
  vtkGetMacro(ScalarFieldId1, int);

  vtkSetMacro(ScalarFieldId2, int);
  vtkGetMacro(ScalarFieldId2, int);

  vtkSetMacro(DistanceType, std::string);
  vtkGetMacro(DistanceType, std::string);

  vtkSetMacro(DistanceFieldName, std::string);
  vtkGetMacro(DistanceFieldName, std::string);

  vtkGetMacro(result, double);

protected:
  ttkLDistance();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string DistanceType{2};
  std::string ScalarField1{};
  std::string ScalarField2{};
  int ScalarFieldId1{0};
  int ScalarFieldId2{1};
  std::string DistanceFieldName{"L2-distance"};
};
