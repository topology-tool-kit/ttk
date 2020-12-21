/// \class ttkPointSetToCurve
/// \ingroup vtk
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
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
#include <ttkPointSetToCurveModule.h>

// VTK includes
#include <ttkAlgorithm.h>

class TTKPOINTSETTOCURVE_EXPORT ttkPointSetToCurve : public ttkAlgorithm {

public:
  static ttkPointSetToCurve *New();
  vtkTypeMacro(ttkPointSetToCurve, ttkAlgorithm);

  vtkSetMacro(CloseCurve, bool);
  vtkGetMacro(CloseCurve, bool);

protected:
  ttkPointSetToCurve();
  ~ttkPointSetToCurve() = default;

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
  bool CloseCurve{false};
};
