/// \class ttkPointsToLine
/// \ingroup vtk
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that reads a Cinema Spec D Database.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \param Output content of the data.csv file of the database in form of a
/// vtkTable

#pragma once

// Module include
#include <ttkPointsToLineModule.h>

// VTK includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

class TTKPOINTSTOLINE_EXPORT ttkPointsToLine : public ttkAlgorithm {

public:
  static ttkPointsToLine *New();
  vtkTypeMacro(ttkPointsToLine, ttkAlgorithm);

  vtkSetMacro(InputOrderingArray, std::string);
  vtkGetMacro(InputOrderingArray, std::string);

protected:
  ttkPointsToLine();
  ~ttkPointsToLine() = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  template <typename VTK_TT>
  void dispatch(std::vector<std::pair<vtkIdType, double>> &storage,
                const VTK_TT *const values,
                const size_t nvalues);

private:
  std::string InputOrderingArray{};
};
