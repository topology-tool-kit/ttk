/// \ingroup vtk
/// \class ttkProjectionFromField
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter which projects a data-set to 2D given two point-data
/// scalar fields to be used as 2D coordinates.
///
/// \param Input Input data-set, with at least two point data scalar fields or
/// texture coordinates (vtkPointSet)
/// \param Output Output projected data-set (vtkPointSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkTextureMapFromField
///
#ifndef _TTK_PROJECTION_FROM_FIELD_H
#define _TTK_PROJECTION_FROM_FIELD_H

// VTK includes
#include <vtkDataArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPointSetAlgorithm.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkProjectionFromField
#else
class ttkProjectionFromField
#endif
  : public vtkPointSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkProjectionFromField *New();

  vtkTypeMacro(ttkProjectionFromField, vtkPointSetAlgorithm);

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

  vtkSetMacro(UComponent, std::string);
  vtkGetMacro(UComponent, std::string);

  vtkSetMacro(VComponent, std::string);
  vtkGetMacro(VComponent, std::string);

  vtkSetMacro(UseTextureCoordinates, bool);
  vtkGetMacro(UseTextureCoordinates, bool);

protected:
  ttkProjectionFromField();

  ~ttkProjectionFromField();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseAllCores;
  ttk::ThreadId ThreadNumber;
  bool UseTextureCoordinates;
  std::string UComponent, VComponent;
  vtkSmartPointer<vtkPoints> pointSet_;

  // base code features
  int doIt(vtkPointSet *input, vtkPointSet *output);

  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};

#endif // _TTK_PROJECTION_FROM_FIELD_H
