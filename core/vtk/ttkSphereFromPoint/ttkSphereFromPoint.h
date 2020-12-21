/// \ingroup vtk
/// \class ttkSphereFromPoint
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter that produces sphere-only glyphs.
///
/// \param Input Input point cloud (vtkDataSet)
/// \param Output Output spheres (vtkPolyData)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
#ifndef _TTK_SPHERE_FROM_POINT_H
#define _TTK_SPHERE_FROM_POINT_H

// VTK includes

// VTK Module
#include <ttkSphereFromPointModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class vtkAppendPolyData;
class vtkSphereSource;

class TTKSPHEREFROMPOINT_EXPORT ttkSphereFromPoint : public ttkAlgorithm {

public:
  static ttkSphereFromPoint *New();

  // macros
  vtkTypeMacro(ttkSphereFromPoint, ttkAlgorithm);

  vtkSetMacro(EndPhi, int);

  vtkSetMacro(EndTheta, int);

  vtkSetMacro(PhiResolution, int);

  vtkSetMacro(Radius, double);

  vtkSetMacro(StartPhi, int);

  vtkSetMacro(StartTheta, int);

  vtkSetMacro(ThetaResolution, int);

  //   /// Over-ride the input data type to vtkDataSet.
  //   int FillOutputPortInformation(int port, vtkInformation *info) override {
  //     info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  //     return 1;
  //   }

protected:
  ttkSphereFromPoint();

  ~ttkSphereFromPoint() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int ThetaResolution{20}, StartTheta{0}, EndTheta{360}, PhiResolution{20},
    StartPhi{0}, EndPhi{180};
  double Radius{0.5};

  vtkAppendPolyData *masterAppender_;
  std::vector<vtkAppendPolyData *> appenderList_;
  std::vector<vtkSphereSource *> sphereList_;
  std::vector<std::vector<vtkDataArray *>> dataArrayList_;
};

#endif // _TTK_SPHERE_FROM_POINT_H
