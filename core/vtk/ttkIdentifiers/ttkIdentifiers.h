/// \ingroup vtk
/// \class ttkIdentifiers
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2015.
///
/// \brief TTK VTK-filter that computes the global identifiers for each vertex
/// and each cell as point data and cell data scalar fields.
///
/// This filter is useful to retrieve the global identifiers of vertices or
/// cells in subsequent filters throughout the VTK pipeline.
///
/// \param Input Input data-set (vtkDataSet)
/// \param Output Output data-set with identifier fields (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
#ifndef _TTK_IDENTIFIERS_H
#define _TTK_IDENTIFIERS_H

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkFiltersCoreModule.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkIdentifiers
#else
class ttkIdentifiers
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkIdentifiers *New();

  vtkTypeMacro(ttkIdentifiers, vtkDataSetAlgorithm);

  vtkSetMacro(CellFieldName, std::string);
  vtkGetMacro(CellFieldName, std::string);

  vtkSetMacro(VertexFieldName, std::string);
  vtkGetMacro(VertexFieldName, std::string);

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

protected:
  ttkIdentifiers();

  ~ttkIdentifiers();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseAllCores;
  ttk::ThreadId ThreadNumber;
  std::string CellFieldName, VertexFieldName;

  // base code features
  int doIt(vtkDataSet *input, vtkDataSet *output);

  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};

#endif // _TTK_IDENTIFIERS_H
