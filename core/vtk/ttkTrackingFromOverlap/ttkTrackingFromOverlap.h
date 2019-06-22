/// \ingroup vtk
/// \class ttkTrackingFromOverlap
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that computes the overlap between labeled
/// vtkPointSets.
///
/// VTK wrapping code for the @TrackingFromOverlap package.
///
/// This filter identifies and tracks labled vtkPointSets across time (and
/// optionally levels) based on spatial overlap, where two points overlap iff
/// their corresponding coordinates are equal. This filter can be executed
/// iteratively and can generate nested tracking graphs.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication: \n
/// 'Nested Tracking Graphs'.
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///
/// \param Input A \b vtkMultiBlockDataSet that holds the labeled \b
/// vtkPointSets and has one of the following forms:\n{t_0,...,t_n} or
/// {l_0:{t_0,...,t_n}, ... , l_m:{t_0,...,t_n}} where \b t_i is the \b
/// vtkPointSet of timestep \b i, and \b l_j is a \b vtkMultiBlockDataSet that
/// holds all timesteps of level \b j. \param Output A \b vtkUnstructuredGrid
/// that represents the (nested) tracking graph embedded in the spatial domain.
///
/// sa ttk::TrackingFromOverlap

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridAlgorithm.h>

// TTK includes
#include <TrackingFromOverlap.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTrackingFromOverlap
#else
class ttkTrackingFromOverlap
#endif
  : public vtkUnstructuredGridAlgorithm,
    public ttk::Wrapper {

public:
  static ttkTrackingFromOverlap *New();
  vtkTypeMacro(ttkTrackingFromOverlap, vtkUnstructuredGridAlgorithm)

    vtkSetMacro(LabelFieldName, string);
  vtkGetMacro(LabelFieldName, string);

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
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkTrackingFromOverlap() {
    SetLabelFieldName("RegionId");

    UseAllCores = false;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkTrackingFromOverlap(){};

  bool UseAllCores;
  int ThreadNumber;

  int reset();

  int packInputData(vtkDataObject *inputDataObject,
                    vtkMultiBlockDataSet *packedData) const;
  int packStreamedData(vtkMultiBlockDataSet *streamedData,
                       vtkMultiBlockDataSet *packedData) const;
  int checkData(vtkMultiBlockDataSet *data);

  int storeStreamedData(vtkMultiBlockDataSet *data);
  int computeNodes(vtkMultiBlockDataSet *data);
  int computeTrackingGraphs(vtkMultiBlockDataSet *data);
  int computeNestingTrees(vtkMultiBlockDataSet *data);
  int computeBranches();

  int meshNestedTrackingGraph(vtkDataObject *trackingGraph);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int LabelDataType;
  string LabelFieldName;
  ttk::TrackingFromOverlap trackingFromOverlap;

  vtkSmartPointer<vtkMultiBlockDataSet> previousIterationData;

  // Containers for nodes and edges
  vector<vector<Nodes>> levelTimeNodesMap; // N
  vector<vector<Edges>> levelTimeEdgesTMap; // E_T
  vector<vector<Edges>> timeLevelEdgesNMap; // E_N

  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};
