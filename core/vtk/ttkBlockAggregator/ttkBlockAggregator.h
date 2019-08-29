/// \ingroup vtk
/// \class ttkBlockAggregator
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2018
///
/// \brief TTK VTK-filter that iteratively appends its input to a
/// vtkMultiBlockDataSet.
///
/// This filter iteratively appends its input as a block to a
/// vtkMultiBlockDataSet.
///
/// \param Input vtkDataObject that will be added as a block (vtkDataObject).
/// \param Output vtkMultiBlockDataSet containing all added blocks
/// (vtkMultiBlockDataSet).

#pragma once

// VTK includes
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkSmartPointer.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkBlockAggregator
#else
class ttkBlockAggregator
#endif
  : public vtkMultiBlockDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkBlockAggregator *New();
  vtkTypeMacro(ttkBlockAggregator, vtkMultiBlockDataSetAlgorithm)

    vtkSetMacro(ForceReset, bool);
  vtkGetMacro(ForceReset, bool);

  vtkSetMacro(FlattenInput, bool);
  vtkGetMacro(FlattenInput, bool);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);
  void SetThreads() {
    threadNumber_
      = !UseAllCores ? ThreadNumber : ttk::OsCall::getNumberOfCores();
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

  int FillInputPortInformation(int port, vtkInformation *info) override {
    if(port < 0 || port > 4)
      return 0;
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
    if(port > 0)
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkBlockAggregator() {
    SetForceReset(false);
    SetFlattenInput(true);

    UseAllCores = false;

    SetNumberOfInputPorts(5);
    SetNumberOfOutputPorts(1);
  }
  ~ttkBlockAggregator(){};

  bool UseAllCores;
  int ThreadNumber;

  int AggregateBlock(vtkDataObject *dataObject, bool useShallowCopy);
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceReset;
  bool FlattenInput;
  vtkSmartPointer<vtkMultiBlockDataSet> AggregatedMultiBlockDataSet;

  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};
