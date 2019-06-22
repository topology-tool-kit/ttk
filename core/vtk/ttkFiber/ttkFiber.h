/// \ingroup vtk
/// \class ttkFiber
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date May 2016.
///
/// \brief TTK VTK-filter for fiber computation on bivariate volumetric data.
///
/// Given a point in the range, this filter computes its fiber (i.e. pre-image)
/// on bivariate volumetric data. The bivariate input data must be provided as
/// two independent scalar fields attached as point data to the input geometry.
///
/// \param Input Input bivariate volumetric data-set, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output Fiber (vtkPolyData)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Fast and Exact Fiber Surface Extraction for Tetrahedral Meshes" \n
/// Pavol Klacansky, Julien Tierny, Hamish Carr, Zhao Geng \n
/// IEEE Transactions on Visualization and Computer Graphics, 2016.
///
/// \sa ttkFiberSurface
///
#ifndef _TTK_FIBER_H
#define _TTK_FIBER_H

// VTK includes -- to adapt
#include <vtkContourFilter.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <Wrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFiber
#else
class ttkFiber
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkFiber *New();

  vtkTypeMacro(ttkFiber, vtkDataSetAlgorithm);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
  }

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  vtkGetMacro(Ucomponent, std::string);
  vtkSetMacro(Ucomponent, std::string);

  vtkGetMacro(Uvalue, double);
  vtkSetMacro(Uvalue, double);

  vtkGetMacro(Vcomponent, std::string);
  vtkSetMacro(Vcomponent, std::string);

  vtkGetMacro(Vvalue, double);
  vtkSetMacro(Vvalue, double);

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }

protected:
  ttkFiber();

  ~ttkFiber();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseAllCores;
  int ThreadNumber;

  double Uvalue, Vvalue;
  std::string Ucomponent, Vcomponent;

  // base code features
  int doIt(vtkDataSet *input, vtkPolyData *output);

  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};

#endif // _TTK_RANGEPOLYGON_H
