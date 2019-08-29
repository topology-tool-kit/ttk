/// \ingroup vtk
/// \class ttkComponentSize
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2015.
///
/// \brief TTK VTK-filter that computes the connected components of a data-set,
/// and that computes their size (number of vertices, number of cells, etc.).
///
/// This filter computes the connected component of a point-set data-set and
/// computes their size (number of vertices, number of cells, etc). The size
/// information is attached on the output, either as point or cell data.
/// The identifier of each connected component is also attached to the geometry
/// as point or cell data.
///
/// This filter is useful when used in conjunction with some thresholding, to
/// only display the largest connected components of a data-set.
///
/// \param Input Input data-set (vtkPointSet)
/// \param Output Output data-set (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
#ifndef _TTK_PROJECTION_FROM_FIELD_H
#define _TTK_PROJECTION_FROM_FIELD_H

// VTK includes
#include <vtkCellData.h>
#include <vtkConnectivityFilter.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPointSetAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// ttk code includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkComponentSize
#else
class ttkComponentSize
#endif
  : public vtkPointSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkComponentSize *New();

  vtkTypeMacro(ttkComponentSize, vtkPointSetAlgorithm);

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

  /// Over-ride the input data type to vtkDataSet.
  int FillOutputPortInformation(int port, vtkInformation *info) override {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }

protected:
  ttkComponentSize();

  ~ttkComponentSize();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseAllCores;
  int ThreadNumber;

  vtkSmartPointer<vtkDoubleArray> vertexNumbers_, cellNumbers_;
  vtkSmartPointer<vtkConnectivityFilter> connectivityFilter_;

  // base code features
  int doIt(vtkPointSet *input, vtkUnstructuredGrid *output);

  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};

#endif // _TTK_PROJECTION_FROM_FIELD_H
