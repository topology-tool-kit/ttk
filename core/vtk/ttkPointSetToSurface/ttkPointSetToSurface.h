/// \class ttkPointSetToSurface
/// \ingroup vtk
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date April 2022
///
/// \brief TTK VTK-filter that reads a Cinema Spec D Database.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \param Output content of the data.csv file of the database in form of a
/// vtkTable
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// Module include
#include <ttkPointSetToSurfaceModule.h>

// VTK includes
#include <ttkAlgorithm.h>

class TTKPOINTSETTOSURFACE_EXPORT ttkPointSetToSurface : public ttkAlgorithm {

public:
  static ttkPointSetToSurface *New();
  vtkTypeMacro(ttkPointSetToSurface, ttkAlgorithm);

protected:
  ttkPointSetToSurface();
  ~ttkPointSetToSurface() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  template <typename VTK_T1, typename VTK_T2>
  void dispatch(std::vector<std::tuple<vtkIdType, double, double>> &storage,
                const VTK_T1 *const values,
                const VTK_T2 *const values2,
                const size_t nvalues);
};
