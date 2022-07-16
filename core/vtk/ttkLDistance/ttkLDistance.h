/// \ingroup vtk
/// \class ttkLDistance
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 26/02/2017
///
/// \brief TTK VTK-filter that wraps the lDistance processing package.
///
/// VTK wrapping code for the ttk::LDistance package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
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

  vtkSetMacro(DistanceType, const std::string &);
  vtkGetMacro(DistanceType, std::string);

  vtkSetMacro(DistanceFieldName, const std::string &);
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
  std::string DistanceType{"2"};
  std::string DistanceFieldName{"L2-distance"};
};
