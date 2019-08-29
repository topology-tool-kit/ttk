/// \ingroup vtk
/// \class ttkGeometrySmoother
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter for geometry smoothing.
///
/// This filter is a dummy example for the development of TTK packages. It
/// smooths an input mesh by average the vertex locations on the link of each
/// vertex.
///
/// \param Input Input mesh (vtkDataSet)
/// \param Output Output mesh (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkScalarFieldSmoother
/// \sa ttk::ScalarFieldSmoother

#ifndef _TTK_GEOMETRY_SMOOTHER_H
#define _TTK_GEOMETRY_SMOOTHER_H

// VTK includes
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <ScalarFieldSmoother.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkGeometrySmoother
#else
class ttkGeometrySmoother
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkGeometrySmoother *New();

  vtkTypeMacro(ttkGeometrySmoother, vtkDataSetAlgorithm);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

  vtkSetMacro(MaskIdentifier, int);
  vtkGetMacro(MaskIdentifier, int);

  vtkSetMacro(ForceInputMaskScalarField, bool);
  vtkGetMacro(ForceInputMaskScalarField, bool);

  vtkSetMacro(InputMask, std::string);
  vtkGetMacro(InputMask, std::string);

protected:
  ttkGeometrySmoother();

  ~ttkGeometrySmoother();

  TTK_SETUP();

private:
  int NumberOfIterations;
  int MaskIdentifier;
  bool ForceInputMaskScalarField;
  std::string InputMask;

  ttk::ScalarFieldSmoother smoother_;
};

#endif // _TTK_GEOMETRY_SMOOTHER_H
