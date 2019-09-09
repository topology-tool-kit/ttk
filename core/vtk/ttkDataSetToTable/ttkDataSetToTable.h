/// \ingroup vtk
/// \class ttkDataSetToTable
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date September 2018
///
/// \brief TTK VTK-filter that creates a vtkTable from a vtkDataSet.
///
/// \param Input vtkDataSet and scalar fields
/// \param Output vtkTable and fields as columns
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
#pragma once

// VTK includes
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

// ttk code includes
#include <Wrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDataSetToTable
#else
class ttkDataSetToTable
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  enum AssociationType { Point = 0, Cell };

  static ttkDataSetToTable *New();
  vtkTypeMacro(ttkDataSetToTable, vtkDataSetAlgorithm)

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

  vtkSetMacro(DataAssociation, int);
  vtkGetMacro(DataAssociation, int);

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
        break;
      default:
        break;
    }

    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {

    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
        break;
      default:
        break;
    }

    return 1;
  }

protected:
  ttkDataSetToTable() {
    UseAllCores = true;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  ~ttkDataSetToTable(){};

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseAllCores;
  int ThreadNumber;
  int DataAssociation;

  int doIt(vtkDataSet *input, vtkTable *output);
  bool needsToAbort() override;
  int updateProgress(const float &progress) override;
};
