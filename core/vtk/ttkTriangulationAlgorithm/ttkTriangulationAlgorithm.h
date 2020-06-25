/// \defgroup vtk vtk
/// \brief The Topology ToolKit - VTK wrapping code for the processing
/// packages.
/// @{
/// \ingroup vtk
/// \class ttkWrapper
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2016.
///

#ifndef _TTK_TRIANGULATION_ALGORITHM_H
#define _TTK_TRIANGULATION_ALGORITHM_H

#include <Wrapper.h>
#include <ttkTriangulation.h>

#include <vtkDataSetAlgorithm.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <macro.h>

#include <ttkTriangulationAlgorithmModule.h>

class TTKTRIANGULATIONALGORITHM_EXPORT ttkTriangulationAlgorithm
  : public vtkDataSetAlgorithm,
    protected ttk::Wrapper {

public:
  static ttkTriangulationAlgorithm *New();

  vtkTypeMacro(ttkTriangulationAlgorithm, vtkDataSetAlgorithm);

protected:
  ttkTriangulationAlgorithm();
  ~ttkTriangulationAlgorithm() override{};

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int RequestDataObject(vtkInformation *request,
                        vtkInformationVector **inputVector,
                        vtkInformationVector *outputVector) override;

private:
  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};

/// @}

#endif /* end of include guard: _TTK_TRIANGULATION_ALGORITHM_H */
